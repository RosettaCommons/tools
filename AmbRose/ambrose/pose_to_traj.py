'''Functions involved in the creation of AMBER simulations based on poses.
Notably, many functions, such as `run_md`_, deal with pure AMBER. `run_md`_ in
particular is, however, the only way to actually run minimizations and
simulations, so it belongs in this family.'''

import re
import tempfile
import os
import subprocess
import sys

# pragma pylint: disable=import-error
import pyrosetta as pr
# pragma pylint: enable=import-error

from . import consts
from . import errors
from . import utils
from . import enums

### Backend utilities

def _heavy_atom_mask(pose):
    mask = pr.rosetta.core.id.AtomID_Map_bool()
    mask.default_value(True)
    size = pose.size()
    mask.resize(size)
    for i in range(1, size+1):
        r = pose.residue(i)
        natoms = r.natoms()
        mask.resize(i, natoms)
        nheavyatoms = r.nheavyatoms()
        for j in range(1, natoms+1):
            if j > nheavyatoms:
                mask.set(pr.AtomID(j, i), False)
    return mask

class _NullContextManager:
    '''Does absolutely nothing when used as a context manager. Equivalent to
    Python 3.7's ``nullcontext`` class from the module ``contextlib``.'''

    def __enter__(self):
        return None
    def __exit__(self, exc_type, exc_value, traceback):
        return False

### Front-facing utilities

def amber_bin():
    '''Get the path to the folder where AMBRose looks for AMBER binaries.'''

    try:
        return os.path.join(os.environ['AMBERHOME'], 'bin')
    except KeyError:
        print('$AMBERHOME is not set in your shell environment! Better hope '
              'all the AMBER executables are in your $PATH. (They should be in '
              '$AMBERHOME/bin.)', file=sys.stderr)
        return ''

def dump_amber_pdb(pose, pdb_path):
    '''From a Pose, dump a PDB suitable for input into AMBER. Currently just
    dumps it normally, but without hydrogens.

    Parameters
    ----------
    pose : rosetta.core.Pose
        The pose to dump.
    pdb_path : str
        The path at which to dump it.
    '''

    # map residue base names to their 3-letter codes in the AMBER naming scheme:
    amber_names = {
        'CYS:disulfide': 'CYX',
        'CYZ':           'CYM',
        'HIS_D':         'HID',
        'HIS':           'HIE',
        'HPR':           'HYP',
    }
    residue_types = {name: pose.conformation() \
                               .residue_type_set_for_conf() \
                               .name_map(name) \
                     for name in amber_names}
    # save original residue names:
    original_names = {name: r.name3() for name, r in residue_types.items()}
    # change the names of all the special residues to the AMBER ones:
    for name, residue_type in residue_types.items():
        residue_type.name3(amber_names[name])
    # actually dump file:
    #pylint: disable=protected-access
    fstream = pr.rosetta.std.ofstream(pdb_path,
                                      pr.rosetta.std._Ios_Openmode._S_out)
    pr.rosetta.core.io.pdb.dump_pdb(pose, fstream, _heavy_atom_mask(pose))
    # change them back:
    for name, residue_type in residue_types.items():
        residue_type.name3(original_names[name])

def dict_to_namelist_str(d, name='cntrl'):
    '''Dumps a single-group FORTRAN 77 NAMELIST for a dict, as a string. The
    string has a trailing newline, to aid concatenation of groups, and to
    make it easier to save as a file. It doesn't do arrays.

    Note that if you just dump the string directly to a file, then it can't
    be used as an mdin NAMELIST file for sander, because sander treats the
    first line as a title, while here the first line is already the group
    name. You need to concatenate a newline to the front of the string, if
    you are going to dump it to a file. The same applies to pmemd.

    Parameters
    ----------
    d : dict
        The dict containing the record names and values for the NAMELIST
        records. Keys should be strings containing valid record names
        (satisfying the regular expression /[A-Za-z][A-Za-z0-9_]{0,30}/), while
        values should be anything that can be cast to a FORTRAN 77 data type that
        isn't an array (ints, floats, strings, booleans).
    name : str, optional
        The name of the group to be output. Must be a valid group name
        (satisfying the regular expression /[A-Za-z][A-Za-z0-9_]{0,30}/).

    Returns
    -------
    str
        A string of your created NAMELIST group. Has a trailing, but not a
        leading newline.
    '''

    valid_id_matcher = re.compile(r'[A-Za-z][A-Za-z0-9_]{0,30}')
    assert valid_id_matcher.fullmatch(name), \
           'Invalid identifier for name of group.'
    assert None not in [valid_id_matcher.fullmatch(key) for key in d.keys()], \
           'Invalid identifier for one or more keys in input dict.'
    to_output = [f'&{name}']
    keys_and_values = tuple(d.items())
    for key, value in keys_and_values[:-1]:
        to_output.append(f'  {key} = {fortran_str(value)},')
    # make sure the last entry doesn't have a comma:
    key, value = keys_and_values[-1]
    to_output.append(f'  {key} = {fortran_str(value)}')
    to_output.extend(['&end', ''])
    return '\n'.join(to_output)

def dict_to_namelist_file(d, path, name='cntrl'):
    '''Creates a single-group FORTRAN 77 NAMELIST file for a given dict of
    key-value pairs ``d`` and a given group name ``name``. The file is saved at
    ``path``. Depends on the function `dict_to_namelist_str()`.

    Parameters
    ----------
    d : dict
        The dict containing the record names and values for the NAMELIST
        records. Keys should be strings containing valid record names
        (satisfying the regular expression /[A-Za-z][A-Za-z0-9_]{0,30}/), while
        values should be anything that can be cast to a FORTRAN 77 data type that
        isn't an array (ints, floats, strings, booleans).
    path : str
        The path at which the output file is to be saved.
    name : str, optional
        The name of the group to be output. Must be a valid group name
        (satisfying the regular expression /[A-Za-z][A-Za-z0-9_]{0,30}/).
    '''

    open(path, 'w').write('\n'+dict_to_namelist_str(d, name))

def fortran_str(x):
    '''Converts ints, floats, bools, and strings into strings of analogous
    FORTRAN 77 literals. Strings must naturally only consist of ASCII, and as a
    further constraint must also only contain printable characters.'''

    if isinstance(x, (int, float)):
        return str(x)
    if isinstance(x, bool):
        if x:
            return '.TRUE.'
        return '.FALSE.'
    if isinstance(x, str):
        to_output = ["'"]
        for c in x:
            if c == "'":
                to_output.extend(("'", "'"))
            else:
                to_output.append(c)
        to_output.append("'")
        to_output = ''.join(to_output)
        assert to_output.isascii() and to_output.isprintable(), \
               'Strings must consist of printable ASCII characters only.'
        return to_output
    raise ValueError('Can only convert ints, floats, bools, and certain '
                     'strings.')

### The money functions

def pose_to_amber_params(pose, crd_path, top_path, log_path='leap.log', *,
                         leap_script_dump_path=None, solvent=None,
                         solvent_shape=enums.SolventShape.OCT, add_ions=True):
    '''Writes a .rst7 input coordinates file and a .parm7 topology file for
    input into sander or pmemd, based on a Rosetta Pose, by invoking LEaP.

    Parameters
    ----------
    pose : rosetta.core.Pose
        The input pose.
    crd_path : str
        Path to the output .rst7 file.
    top_path : str
        Path to the output .parm7 file.
    log_path : str, optional
        Path to the output LEaP log file. Default is ``leap.log`` in the current
        directory. Due to a how LEaP's header scripts for proteins function,
        this command actually barely does anything at all, at the moment; most
        of the log will still be dumped in the current working directory under
        the filename ``leap.log``, regardless of your settings.
    leap_script_dump_path : str or None, optional
        Path to dump the generated LEaP script at. Default is None.
    solvent : Solvent or None, optional
        What explicit solvent to use, if any. See the `Solvent`_ class for your
        options. None by default.
    solvent_shape : SolventShape, optional
        What shape to use for the explicit solvent box, if any. See the
        `SolventShape`_ class for your options. ``OCT`` by default.
    add_ions : bool, optional
        Whether to add neutralizing sodium or chloride ions in the case that
        your simulation is explicit-solvent and the charge of your protein is
        nonzero. True by default.
    '''

    # .exists() is appropriate here rather than .isfile(), because we don't want
    # to remove it even if it's a directory. I don't know what LEaP does with it
    # if it already exists as a directory, however...
    leap_log_exists_in_current_dir = os.path.exists('leap.log')
    if leap_script_dump_path is not None:
        leap_in_fd, leap_in_path = leap_script_dump_path, leap_script_dump_path
    else:
        leap_in_fd, leap_in_path = tempfile.mkstemp()
    try:
        with tempfile.NamedTemporaryFile() as pdb_f:
            ## dump the PDB we're going to feed to LEaP
            dump_amber_pdb(pose, pdb_f.name)

            ## generate and write the LEaP instructions
            leap_command_str = \
                open(os.path.join(consts.AMBROSE_DIR,
                                  'pose-to-traj.in.prototype')).read()
            if solvent:
                solvent_str, load_solvent_str = \
                    pose_to_amber_params.SOLVENTS[solvent]
            disulfide_residues = []
            for i in range(1, pose.size()+1):
                if pose.residue(i).name() == 'CYS:disulfide':
                    disulfide_residues.append(i)
            disulfides = []
            while disulfide_residues:
                # TODO: Rewrite this to not run in O(nlogn) time; it's trash.
                our_r = disulfide_residues.pop()
                other_r = None
                for i in disulfide_residues:
                    if pose.residue(our_r).is_bonded(i):
                        other_r = i
                        break
                else:
                    raise errors.TopologyTypeError(
                        'Cystine is bonded to something other than a cystine. '
                        'AMBRose does not support noncanonical residues yet. '
                        'If you think that this is the only error preventing '
                        'your project from working properly, please raise an '
                        'issue on the rosetta/tools Github and give it the '
                        '"ambrose" tag. I\'ll see what I can do about it.')
                disulfides.append((our_r, other_r))
                disulfide_residues.remove(other_r)
            disulfides_str = '\n'.join(f'bond struct.{i}.SG struct.{j}.SG' \
                                       for i, j in disulfides)
            solvent_shape_str = \
                pose_to_amber_params.SOLVENT_SHAPES[solvent_shape]
            # The following code will replace each instance of each key in the
            # leap command str with its corresponding value. The keys must
            # not be prefixes of each other.
            replacements = {'+LOG-FILE+':
                                log_path,
                            '+LOAD-SOLVENT+':
                                load_solvent_str if solvent else '',
                            '+PDB-PATH+':
                                pdb_f.name,
                            '+DISULFIDES+':
                                disulfides_str,
                            '+SOLVATEP+':
                                '' if solvent else '#',
                            '+SOLVENT-SHAPE+':
                                solvent_shape_str if solvent else '',
                            '+SOLVENT+':
                                solvent_str if solvent else '',
                            '+ADD-IONS-P+':
                                '' if add_ions and solvent else '#',
                            '+OUT-TOP-PATH+':
                                top_path,
                            '+OUT-CRD-PATH+':
                                crd_path}
            matcher = re.compile(
                '(' + '|'.join(map(re.escape, replacements.keys())) + ')')
            leap_command_str = matcher.sub(
                lambda mo: replacements[mo.group(1)],
                leap_command_str)
            open(leap_in_fd, 'w').write(leap_command_str)

            ## get LEaP to write the file
            subprocess.run([os.path.join(amber_bin(), 'tleap'),
                            '-f', leap_in_path])

    finally:
        if leap_script_dump_path is None:
            try:
                os.remove(leap_in_path)
            except FileNotFoundError:
                pass
        if not leap_log_exists_in_current_dir:
            try:
                os.remove('leap.log')
            except FileNotFoundError:
                pass
## the "load" commands for different solvents, used to construct the LEaP
## instructions
pose_to_amber_params.SOLVENTS = \
    {enums.Solvent.SPCE_WATER:
     ('SPCBOX',
      'source leaprc.water.spce'),
     enums.Solvent.TIP3P_WATER:
     ('TIP3PBOX',
      'source leaprc.water.tip3p'),
     enums.Solvent.METHANOL:
     ('MEOHBOX',
      'loadoff solvents.lib\nloadamberparams frcmod.meoh')}
pose_to_amber_params.SOLVENT_SHAPES = \
    {enums.SolventShape.BOX: 'box',
     enums.SolventShape.OCT: 'oct'}

def run_md(executable='pmemd.cuda', *, overwrite=False, mdin=None,
           mdin_dict=None, mdout=None, prmtop=None, inpcrd=None, restrt=None,
           refc=None, mdcrd=None, mdinfo=None):
    '''Runs an AMBER MD program (pmemd.cuda by default) with the given arguments
    with the additional feature that the mdin file, which normally gives the
    simulation parameters, can be replaced with a dict specifying the same
    parameters. All arguments are optional, as in the sander executable; if a
    value is needed but not specified, sander will assume the default (which is
    the argument's name in this case).

    If you don't want to make the mdin_dict parameters by hand, the templates_
    submodule can synthesize some common formulas for you.

    Parameters
    ----------
    executable : str, optional
        The AMBER executable to run the MD with. pmemd.cuda by default.
    overwrite : bool, optional
        Whether to overwrite the output files if already present. Equivalent to
        the ``-O`` flag in the executable.
    mdin : str or None, optional
        The path to the FORTRAN NAMELIST file that contains the simulation
        parameters. Overridden if parameter mdin_dict is not None. Equivalent to
        the ``-i`` flag in the executable.
    mdin_dict : dict or None, optional
        A dict encoding the FORTRAN NAMELIST file that would contain the
        simulation parameters. The group name is always ``cntrl`` (hence Ewald
        is currently not supported). The keys and values represent data fields
        and their values respectively. If both this and ``mdin`` is set, then
        the corresponding ``mdin`` file generated from the ``mdin_dict`` will be
        stored at the path given for ``mdin``.
        See the `templates` submodule for functions that automatically
        produce ``mdin_dict``s for different situations.
    mdout : str or None, optional
        The path to the log file where each ``mdinfo`` is recorded. Equivalent
        to the ``-o`` flag in the executable.
    prmtop : str or None, optional
        The path to the topology file of the structures involved in the
        simulation. Equivalent to the ``-p`` flag in the executable.
    inpcrd : str or None, optional
        The path to the file containing the initial coordinates for the
        simulation. Equivalent to the ``-c`` flag in the executable.
    restrt : str or None, optional
        The path to the file containing the most recently computed coordinates
        and velocities for the simulation. Equivalent to the ``-r`` flag in the
        executable.
    refc : str or None, optional
        The path to the file containing the reference coordinates for coordinate
        constraints, if any. Equivalent to the ``-ref`` flag in the executable.
    mdcrd : str or None, optional
        The path to the file where the coordinates are output during the
        simulation. Equivalent to the ``-x`` flag in the executable.
    mdinfo : str or None, optional
        The path to the file where some high-level data about how the simulation
        is going is saved every so often. Equivalent to the ``-inf`` flag in the
        executable.

    Returns
    -------
    int
        Return value of the AMBER executable run.
    '''

    with tempfile.NamedTemporaryFile() if mdin is None \
                                          and mdin_dict is not None \
         else _NullContextManager() as mdin_tempfile:
        if mdin_dict is not None:
            if mdin is None:
                mdin = mdin_tempfile.name
            dict_to_namelist_file(mdin_dict, mdin, name='cntrl')
        flags_to_variables = {
            '-i' : mdin,
            '-o' : mdout,
            '-p' : prmtop,
            '-c' : inpcrd,
            '-r' : restrt,
            '-ref' : refc,
            '-x' : mdcrd,
            '-inf' : mdinfo
        }
        arguments = ['-O'] if overwrite else []
        for key, value in flags_to_variables.items():
            if value is not None:
                arguments.extend((key, value))
        command = [os.path.join(amber_bin(), executable), *arguments]
        utils.debug_print(f'About to run command: {command}')
        result = subprocess.run(command)
        utils.debug_print(f'{executable} run with result: {result}')
        return result
