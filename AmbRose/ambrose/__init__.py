'''Gonna be part of Ambrose. Attempts to import a pytraj Trajectory as a
pyrosetta Pose.'''

# pragma pylint: disable=bad-whitespace

import re
import datetime
import collections
import weakref
import copy
import enum
import itertools
import tempfile
import os
import subprocess
import sys
import warnings

from functools import reduce

# pragma pylint: disable=import-error
import pytraj as pt
import pyrosetta as pr
from pyrosetta.rosetta.core.io import AtomInformation
# pragma pylint: enable=import-error

import pose_selectors

AMBROSE_DIR = os.path.dirname(__file__)
DEBUG = False

### Helper functions and classes for dealing with Rosetta and AMBER
## Get certain Rosetta singleton objects

def _chemical_manager():
    return pr.rosetta.core.chemical.ChemicalManager.get_instance()

def _global_element_set():
    if not hasattr(_global_element_set, 'es'):
        _global_element_set.es = _chemical_manager().element_set('default')
    return _global_element_set.es

def _global_residue_type_set():
    if not hasattr(_global_residue_type_set, 'rts'):
        _global_residue_type_set.rts = \
            _chemical_manager().residue_type_set('fa_standard')
    return _global_residue_type_set.rts

## Work with Rosetta ResidueType objects

def _get_pr_res(name):
    return _global_residue_type_set().get_representative_type_base_name(name)

def _patch_pr_res(res, patch_name):
    return _global_residue_type_set().patch_map()[patch_name][1].apply(res)

## Work with Rosetta Pose objects

def _check_pose_convertibility(pose, crd_path, top_path):
    '''Checks whether a Pose changes its topology upon being converted to AMBER
    coordinates and back, given the pose and the paths to the coordinates and
    topology files it's been converted to.

    Parameters
    ----------
    pose : rosetta.core.Pose
        Pose to check.
    crd_path : str
        Path to coordinates file Pose has been converted to.
    top_path : str
        Path to topology file Pose has been converted to.
    '''

    ambered_pose = TrajToPoses(pt.iterload(crd_path, top_path))[0]
    # pylint: disable=no-member
    if ambered_pose.size() > pose.size():
        bad_residue = None
        for i in range(1, ambered_pose.size()+1):
            if ambered_pose.residue(i).type() != pose.residue(i).type():
                bad_residue = f'{i} {ambered_pose.residue(i).name()}'
                break
        raise TopologySizeError('Pose changed size when piped through '
                                'AMBER! It somehow got bigger, which '
                                'should be impossible. Please open an '
                                'issue on the Github repo and describe '
                                'your error there in detail. Here\'s the '
                                '(first) residue that sprouted out of '
                                'thin air:\n  ' + bad_residue)
    if ambered_pose.size() < pose.size():
        bad_residue = None
        for i in range(1, pose.size()+1):
            if ambered_pose.residue(i).type() != pose.residue(i).type():
                bad_residue = f'{i} {pose.residue(i).name()}'
                break
        raise TopologySizeError('Pose changed size when piped through '
                                'AMBER! Here\'s the (first) residue that '
                                'went missing:\n  ' + bad_residue)
    for i in range(1, pose.size()+1):
        bad_residue = None
        if ambered_pose.residue(i).type() != pose.residue(i).type():
            bad_residue = f'{i} {pose.residue(i).name()}'
            break
    if bad_residue is not None:
        raise TopologyTypeError('Pose changed at least one residue\'s type '
                                'when piped through AMBER! Here\'s the '
                                'offending residue:\n  ' + bad_residue)
    # pylint: enable=no-member


def _transfer_xyz_from_pose_to_pose(source_pose, dest_pose):
    '''Copy coordinates from one pose to another, assuming that they have the
    exact same topology (we don't bother to check).'''

    for i in range(1, dest_pose.size()+1):
        dest_r = dest_pose.residue(i)
        source_r = source_pose.residue(i) # pylint: disable=no-member
        for j in range(1, dest_r.natoms()):
            dest_r.set_xyz(j, source_r.xyz(j))

## Fake Mover class

## Misc Rosetta stuff

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

## Misc AMBER stuff

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

    fstream = pr.rosetta.std.ofstream(pdb_path,
                                      pr.rosetta.std._Ios_Openmode._S_out)
    pr.rosetta.core.io.pdb.dump_pdb(pose, fstream, _heavy_atom_mask(pose))

def dict_to_namelist_str(d, name='cntrl'):
    '''Dumps a single-group FORTRAN 77 NAMELIST for a dict, as a string. The
    string has a trailing newline, to aid concatenation of groups, and to
    make it easier to save as a file. It doesn't do arrays.

    Note that if you just dump the string directly to a file, then it can't
    be used as an mdin NAMELIST file for sander, because sander treats the
    first line as a title, while here the first line is already the group
    name. You need to concatenate a newline to the front of the string, if
    you are going to dump it to a file.

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
    to_output.extend(['&', ''])
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

## Misc other stuff

def _debug_print(value, *more_values, sep=' ', end='\n', file=sys.stdout,
                 flush=True):
    if DEBUG:
        print(value, *more_values, sep=sep, end=end, file=file, flush=flush)

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
                x.extend(("'", "'"))
            else:
                x.append(c)
        to_output.append("'")
        to_output = ''.join(to_output)
        assert to_output.isascii() and to_output.isprintable(), \
               'Strings must consist of printable ASCII characters only.'
        return to_output
    raise ValueError('Can only convert ints, floats, bools, and certain '
                     'strings.')

def _timestamp():
    '''Gets timestamp in standardized format.'''

    return datetime.datetime.now().isoformat()

### Supporting classes

## Misc

class _NullContextManager:
    '''Does absolutely nothing when used as a context manager. Equivalent to
    Python 3.7's ``nullcontext`` class from the module ``contextlib``.'''

    def __enter__(self):
        return None
    def __exit__(self, exc_type, exc_value, traceback):
        return False

## Classes supporting TrajToPoses

class _AtomRecord:
    '''Describes the key properties of an atom in a PyRosetta ResidueType used
    to construct AtomInformation records, and its mapping onto an absolute
    index in a Topology.

    Parameters
    ----------
    parent_res : :obj:`_ResidueMap` or None
        The :obj:`_ResidueMap` that this ``_AtomRecord`` belongs to. Can be None,
        but then you can't call `atom_information`.
    name : str
        The 4-character name of the atom in its Rosetta ``ResidueType``.
    element : str
        The short string representing the element of this atom in Rosetta,
        usually equivalent to the atom's chemical symbol. (See
        ``core::chemical::Elements`` in Rosetta.)
    basep : bool
        Whether or not this atom appears in the unpatched version of its parent
        residue. Used by some ``set_indices`` methods of :obj:`_ResidueMap`.

    Attributes
    ----------
    parent_res : :obj:`_ResidueMap` or None
        The :obj:`_ResidueMap` that this ``_AtomRecord`` belongs to. Can be None,
        but then you can't call `atom_information`.
    name : str
        The 4-character name of the atom in its Rosetta ``ResidueType``.
    h_number : int or None
        The number that appears at the beginning of this atom's name, if it is a
        hydrogen numbered that way.  Used by some ``set_indices`` methods of
        :obj:`_ResidueMap`.
    coda : str
        The atom's name without its `h_number` and surrounding spaces. Equal to
        `name` with its spaces stripped if `h_number` is None.
    pt_i : int or None
        The absolute index of the atom that this atom corresponds to in a
        certain pytraj trajectory. Must be assigned externally, by convention by
        ``set_indices`` methods of its parent :obj:`_ResidueMap`.
    element : str
        The short string representing the element of this atom in Rosetta,
        usually equivalent to the atom's chemical symbol. (See
        ``core::chemical::Elements`` in Rosetta.)
    basep : bool
        Whether or not this atom appears in the unpatched version of its parent
        residue. Used by some ``set_indices`` methods of :obj:`_ResidueMap`.
    '''

    def __init__(self, parent_res, name, element, basep=True):
        self.parent_res = parent_res
        self.name     = name
        self.h_number = int(name[0]) if name[0] in '123456789' else None
        self.coda     = name[1:].strip(' ') if self.h_number is not None \
                        else name.strip(' ')
        self.pt_i     = None
        self.element  = element
        self.basep    = basep
    @property
    def atom_information(self):
        '''Returns a Rosetta ``AtomInformation`` object based on this record.

        Returns
        -------
        rosetta.core.io.AtomInformation'''

        if self.parent_res is None:
            raise NotImplementedError('parent_res cannot be None if '
                                      'atom_information() is called!')
        ai = AtomInformation()
        ai.resName   = self.parent_res.name3
        ai.name      = self.name
        ai.chainID   = self.parent_res.chain_id
        ai.occupancy = 1.0
        ai.segmentID = '    '
        ai.element   = self.element
        return ai

class _ResidueMap:
    '''Describes the key properties of a ResidueType used to construct
    :obj:`_AtomRecord`s, and a mapping of its atoms onto absolute indices in a
    pytraj ``Topology``. This object can be indexed and iterated over to
    retreive its :obj:`_AtomRecord`s, but not sliced.

    Parameters
    ----------
    name : str
        The base name of the corresponding ``ResidueType`` in Rosetta.
    patches : tuple
        A tuple of names of patches to apply to the base residue to get the
        final corresponding ``ResidueType``.
    first_atom_i : int
        The index of the first atom of the residue that this corresponds to in
        pytraj.
    always_base_p : bool, optional
        Whether to treat all atoms as base atoms, rather than just the ones that
        appear in the base ``ResidueType``. This has significance for certain
        ``set_indices`` methods.

    Attributes
    ----------
    name : str
        The base name of the corresponding ``ResidueType`` in Rosetta.
    patches : tuple
        A tuple of names of patches to apply to the base residue to get the
        final corresponding ``ResidueType``.
    first_atom_i : int
        The index of the first atom of the residue that this corresponds to in
        pytraj.
    chain_id : str
        A single printable ASCII character representing the chain to which this
        residue belongs. Currently always equal to ``@`` (0x40), but it is
        planned to break non-contiguous topologies up into multiple chains
        eventually, at which point chain IDs will count up from 0x40 to 0xFF.
    name3 : str
        The three-letter residue name for this residue in Rosetta.
    atoms : list(:obj:`_AtomRecord`)
        A list of the :obj:`_AtomRecord`s of this residue.
    __names_to_indices : dict
        A dict that maps the 4-character names of the atoms in the corresponding
        Rosetta ``ResidueType`` to their indices in `atoms`.
    '''

    def __init__(self, name, patches, first_atom_i, always_base_p=False):
        self.name = name
        self.patches = patches
        self.first_atom_i = first_atom_i
        self.chain_id = '@'

        base_res_type = _get_pr_res(name)
        res_type      = reduce(_patch_pr_res, patches, base_res_type)

        self.name3 = res_type.name3()
        self.atoms = []
        self.__names_to_indices = {}
        for i in range(1, res_type.natoms()):
            atom_name = res_type.atom_name(i)
            self.__names_to_indices[atom_name] = i
            self.atoms.append(
                _AtomRecord(self,
                            atom_name,
                            str(res_type.atom(i).element())[9:],
                            basep=(always_base_p or \
                                   base_res_type.has(atom_name))))
        self.atoms = tuple(self.atoms)

    def set_indices_by_h_number_scheme(self, pt_names_to_indices):
        '''Set the pt_i of the base atoms, based on a map of pytraj atom names
        to absolute pt_is, via a specific algorithm that works for canonical
        amino acids and canonical nucleic acids (and hydroxyproline, and not
        much else).

        Parameters
        ----------
        pt_names_to_indices : dict
            A dict that maps the names of atoms in a pytraj residue to their
            absolute indices in their trajectory.
        '''

        # collection of codas of residues that have h numbers, mapped to a list
        # of the atom records that have it:
        codas_with_h_numbers = collections.defaultdict(list)
        # write indices of residues that don't have h numbers
        for atom_record in self:
            if atom_record.basep:
                ## Mark atoms with h numbers and don't process them:
                if atom_record.h_number is not None:
                    codas_with_h_numbers[atom_record.coda].append(atom_record)
                    continue

                ## Handle edge cases:
                # rosetta prolines have an extra virtual nitrogen at the
                # same spot as the backbone nitrogen:
                if self.name in ('PRO', 'HPR') and atom_record.name == ' NV ':
                    atom_record.h_number = None
                    atom_record.coda = 'N'

                ## Actually retrieve index:
                try:
                    atom_record.pt_i = pt_names_to_indices[atom_record.coda]
                except KeyError:
                    # indicates that a mapping couldn't be established;
                    # hopefully whatever is processing the _ResidueMap will
                    # be able to deal with the missing mapping
                    atom_record.pt_i = None

        # maps codas to starting h number in pt residue (pr residue h
        # numbers assumed to always start at 1):
        pt_starting_h_numbers = \
            {coda: min(int(name[-1]) \
                       for name in pt_names_to_indices.keys() \
                       if name[:-1] == coda and name[-1].isdigit())
             for coda in codas_with_h_numbers.keys()}
        # write indices of residues that DO have h numbers:
        for coda, atom_records in codas_with_h_numbers.items():
            for atom_record in atom_records:
                atom_record.pt_i = \
                    pt_names_to_indices[
                        coda + \
                        str(atom_record.h_number - 1 + \
                            pt_starting_h_numbers[atom_record.coda])]
    def set_indices_by_name_relative(self, indices_dict):
        '''Provided a dict that maps Rosetta atom names to an offset of its pt_i
        from `first_atom_i` (the index of the first atom in the pytraj residue),
        set the pt_is of each named atom to the appropriate absolute index.

        Parameters
        ----------
        indices_dict : dict
            A dict mapping the Rosetta names of atoms in this residue to the
            indices of their corresponding atoms in a pytraj trajectory,
            relative to `first_atom_i`.
        '''

        for atom_name, i in indices_dict.items():
            self[atom_name].pt_i = self.first_atom_i + i
    def __len__(self):
        return len(self.atoms)
    def __iter__(self):
        for i in range(len(self)):
            yield self[i]
    def __getitem__(self, key):
        if isinstance(key, slice):
            raise NotImplementedError
        if isinstance(key, str):
            return self.atoms[self.__names_to_indices[key]]
        return self.atoms[key]

class _StateID:
    SALT = '___StateID_'
    def __init__(self, init_id=0):
        self.id = init_id
    def __hash__(self):
        return hash(_StateID.SALT + str(self.id))
    def __eq__(self, other):
        if isinstance(other, _StateID):
            return self.id == other.id
        return False
    def next(self):
        return _StateID(self.id + 1)

class _State:
    next_id = _StateID(0)
    def __init__(self, id=None):
        if id is not None:
            self.id = id
        else:
            self.id = _State.next_id
            _State.next_id = _State.next_id.next()
    def __hash__(self):
        return hash(self.id)
    def __eq__(self, other):
        return self.id == other.id

class _Transition:
    def __init__(self, next_state, emitted=()):
        self.next_state = next_state
        self.emitted    = emitted
    def __hash__(self):
        return hash(self.next_state, self.emitted)

class _TopologyParser:
    '''A finite state machine that parses residues from a pytraj topology and
    builds a list of _ResidueMap objects.'''

    # Canonical amino acid residue name sets
    AA = frozenset({'ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'CYM', 'CYS', 'CYX',
                    'GLN', 'GLU', 'GLH', 'GLY', 'HID', 'HIE', 'HIP', 'HYP',
                    'ILE', 'LEU', 'LYN', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                    'THR', 'TRP', 'TYR', 'VAL'})
    NAA  = frozenset('N'+aa for aa in AA)
    CAA  = frozenset('C'+aa for aa in AA)
    RNA  = frozenset({'A', 'C', 'G', 'U'})
    RNA3 = frozenset(base+'3' for base in RNA)
    RNA5 = frozenset(base+'5' for base in RNA)
    DNA  = frozenset({'DA', 'DC', 'DG', 'DT'})
    DNA3 = frozenset(base+'3' for base in DNA)
    DNA5 = frozenset(base+'5' for base in DNA)
    NA   = RNA | DNA
    NA3  = RNA3 | DNA3
    NA5  = RNA5 | DNA5

    # Auxiliary data for map_atoms
    PT_TO_PR_RES_NAMES = {
        'ASH':'ASP',
        'CYM':'CYS',
        'CYX':'CYS',
        'GLH':'GLU',
        'HID':'HIS',
        'HIE':'HIS',
        'HIP':'HIS',
        'HYP':'HPR',
        'LYN':'LYS',
        'A'  :'RAD',
        'A3' :'RAD',
        'A5' :'RAD',
        'C'  :'RCY',
        'C3' :'RCY',
        'C5' :'RCY',
        'G'  :'RGU',
        'G3' :'RGU',
        'G5' :'RGU',
        'U'  :'URA',
        'U3' :'URA',
        'U5' :'URA',
    }
    TOKENS_TO_PATCHES = {
        'NAA':   ('NtermProteinFull',),
        'ACE_AA':('N_acetylated',),
        'CAA':   ('CtermProteinFull',),
        'AA_NME':('C_methylamidated',),
        'NA3':   ('3prime_end_OH',),
        'NA5':   ('5prime_end_OH',), # TODO: find out what patch it actually is
    }

    # States:
    BEGIN_STATE  = _State()
    ACE_STATE    = _State()
    PRE_AA_STATE = _State()
    AA_STATE     = _State()
    NA_STATE     = _State()
    END_STATE    = _State()

    def __init__(self):
        ## build transitions:
        # brevity:
        BEGIN  = _TopologyParser.BEGIN_STATE
        ACE    = _TopologyParser.ACE_STATE
        PRE_AA = _TopologyParser.PRE_AA_STATE
        AA     = _TopologyParser.AA_STATE
        NA     = _TopologyParser.NA_STATE
        END    = _TopologyParser.END_STATE
        # emitted objects are tuples of "residue tokens" that are used to build
        # the actual _ResidueMaps
        self.state = BEGIN
        self.transitions = \
            {BEGIN:  {_TopologyParser.NAA: (PRE_AA, ('NAA',)),
                      _TopologyParser.AA:  (PRE_AA, ('NAA',)),
                      frozenset({'ACE'}):  (ACE,),
                      _TopologyParser.NA:  (NA,     ('NA',)),
                      _TopologyParser.NA3: (NA,     ('NA3',)),
                      _TopologyParser.NA5: (NA,     ('NA5',))},
             ACE:    {_TopologyParser.AA:  (PRE_AA, ('ACE_AA',))},
             PRE_AA: {_TopologyParser.AA:  (AA,)},
             AA:     {_TopologyParser.NAA: (PRE_AA, ('CAA','NAA')),
                      _TopologyParser.AA:  (AA,     ('AA',)),
                      _TopologyParser.CAA: (END,    ('AA', 'CAA')),
                      frozenset({'NME'}):  (END,    ('AA_NME',))},
             NA:     {_TopologyParser.NA:  (NA,     ('NA',)),
                      _TopologyParser.NA3: (NA,     ('NA3',)),
                      _TopologyParser.NA5: (NA,     ('NA5',))},
             END:    {}}
        self.default_transitions = \
            {BEGIN:  (BEGIN, ('SKIP',)),
             ACE:    (BEGIN, ('SKIP',)),
             PRE_AA: (BEGIN, ('SKIP',)),
             AA:     (END,   ('CAA', 'SKIP')),
             NA:     (END,   ('SKIP',)),
             END:    (END,   ('SKIP',))}
    def transit(self, received_id):
        '''Perform a state transition, with received_id as the key that maps to
        the transition.'''

        for id_set, possible_transition in self.transitions[self.state].items():
            if received_id in id_set:
                our_transition = _Transition(*possible_transition)
                break
        else: # if we went through the whole loop without breaking...
            our_transition = _Transition(*self.default_transitions[self.state])
        self.state = our_transition.next_state
        return our_transition.emitted

    def parse_topology(self, top):
        '''Returns a list of `_ResidueMap`_ objects, which contain the
        information needed to build core::chemical::Residue objects from the
        list of atoms in the given topology, namely ResidueType names, any
        Patches, and a map of each ResidueType's atoms to an atom index in the
        given topology.'''

        maps = [] # collect output maps here
        pt_res_i = 0
        pt_res = top.residue(pt_res_i)
        pt_atom_i = pt_res.first_atom_index

        def consume_residue():
            '''Updates pt_res_i, pt_atom_i, and n_res_consumed to the next
            residue.'''
            nonlocal pt_res_i
            nonlocal pt_res
            nonlocal pt_atom_i
            pt_res_i += 1
            if pt_res_i < top.n_residues:
                pt_res = top.residue(pt_res_i)
                pt_atom_i = pt_res.first_atom_index
            _debug_print('residue consumed:', pt_res_i, pt_res)

        for residue in top.residues:
            _debug_print('residue:', residue)
            tokens = list(self.transit(residue.name))
            if pt_res_i == top.n_residues-2 and tokens[-1] == 'AA':
                tokens.append('CAA')
            for token in tokens:
                _debug_print('token:', token)
                if token == 'SKIP':
                    consume_residue()
                    continue

                naa_or_caa_p = token in ('NAA', 'CAA')
                na3_or_na5_p = token in ('NA3', 'NA5')

                ## determine name of new residue, name
                pt_name = pt_res.name
                name = _TopologyParser.PT_TO_PR_RES_NAMES.get(pt_name, pt_name)
                if naa_or_caa_p and len(name) == 4 and name[0] == 'N':
                    name = name[1:]

                ## determine patches to be applied, patches
                patches = _TopologyParser.TOKENS_TO_PATCHES.get(token, ())

                ## construct _ResidueMap records:
                res_map = _ResidueMap(name, patches, pt_atom_i,
                                      always_base_p=(naa_or_caa_p or \
                                                     na3_or_na5_p))
                # assign indices for acetyl atoms
                if token == 'ACE_AA':
                    res_map.set_indices_by_name_relative(
                        {' CP ': 4,
                         ' CQ ': 1,
                         ' OCP': 5})
                    consume_residue()
                # assign indices for atoms of base residue
                if token in ('NAA', 'ACE_AA', 'AA', 'CAA', 'AA_NME', 'NA',
                             'NA3', 'NA5'):
                    res_map.set_indices_by_h_number_scheme(
                        {top[i].name: i \
                         for i \
                         in range(pt_atom_i, pt_atom_i + pt_res.n_atoms)})
                    consume_residue()
                # assign indices for C-terminal N-methylamine atoms
                if token == 'AA_NME':
                    res_map.set_indices_by_name_relative(
                        {' NR ': 0,
                         ' CS ': 2,
                         ' HR ': 1})
                    consume_residue()

                maps.append(res_map)
        return tuple(maps)

## Classes supporting the _AMBERMovers

class _NotAMover:
    '''Mimics the Mover class from Rosetta. Nearly every method raises a
    NotImplementedError; it is therefore more of an interface than a class.
    All getters and setters are implemented as Python properties, with "get"
    and "set" removed from their names, if applicable. The methods that don't
    raise NotImplementedError have docstrings.'''

    def __init__(self):
        raise NotImplementedError
    @property
    def name(self):
        '''Name of Mover.'''
        return 'UNDEFINED NAME'
    @property
    def self_ptr(self):
        '''Pointer to Mover.'''
        return self
    @property
    def self_weak_ptr(self):
        '''Weak pointer to Mover.'''
        return weakref.ref(self)
    @staticmethod
    def register_options():
        raise NotImplementedError
    def apply(self, pose):
        raise NotImplementedError
    def test_move(self, pose):
        '''Copy of `apply`_, with possibility of overriding for testing.'''
        self.apply(pose)
    @property
    def type(self):
        raise NotImplementedError
    @type.setter
    def type(self, value):
        raise NotImplementedError
    @property
    def input_pose(self):
        raise NotImplementedError
    @input_pose.setter
    def input_pose(self, value):
        raise NotImplementedError
    @property
    def native_pose(self):
        raise NotImplementedError
    @native_pose.setter
    def native_pose(self, value):
        raise NotImplementedError
    def clone(self):
        '''Returns copy of Mover.'''
        return copy.copy(self)
    def parse_my_tag(self, tag, data, filters, movers, pose):
        raise NotImplementedError
    def show(self, ostream=sys.stdout):
        '''Prints a representation of the Mover.'''
        print(repr(self), file=ostream)

class _AMBERMover(_NotAMover):
    '''Not a real mover; rather, a superclass for all "movers" using AMBER.
    The apply_ method actually does something: it minimizes the pose for
    500 steps in AMBER. This doesn't change the final state of the pose, but
    it doesn't delete the files, either. This is meant to serve as the first
    step of any AMBER operation, and is ideally called as a super before the
    subclass's apply. See its docstring for details on the kinds of files it
    creates.'''

    LEAP_SUFFIXES = {
        'crd_path' : '-inpcrd.rst7',
        'top_path' : '-prmtop.parm7',
        'log_path' : '-leap.log'
    }
    MIN_SUFFIXES = {
        'prmtop' : '-prmtop.parm7',
        'inpcrd' : '-inpcrd.rst7',
        'restrt' : '-min-restrt.nc',
        'mdinfo' : '-min-mdinfo',
        'mdin'   : '-min-mdin.in',
        'mdout'  : '-min-mdout.out',
        'refc'   : '-refc.rst7'
    }
    MD_SUFFIXES = {
        'prmtop' : '-prmtop.parm7',
        'inpcrd' : '-min-restrt.nc',
        'restrt' : '-restrt.nc',
        'mdcrd'  : '-mdcrd.nc',
        'mdinfo' : '-mdinfo',
        'mdin'   : '-mdin.in',
        'mdout'  : '-mdout.out',
        'refc'   : '-refc.rst7'
    }

    def __init__(self, working_dir='', prefix=None, solvent=None,
                 solvent_shape=SolventShape.OCT, cutoff=15., ref_pose=None,
                 cst_mask=None, cst_weight=None):
        self._working_dir = working_dir
        self._prefix = prefix
        self._ref_pose = ref_pose
        self._solvent = solvent
        self._solvent_shape = solvent_shape
        wet = solvent is not None
        self._amber_executable = 'pmemd' if wet else 'sander'
        common_dict = {'igb': int(not wet),
                       'ntb': int(wet),
                       'cut': cutoff}
        if cst_mask is not None:
            common_dict['restraintmask'] = cst_mask
        if cst_weight is not None:
            common_dict['restraint_wt'] = cst_weight
        self._min_mdin_dict = {'imin': 1, **common_dict}
        self._mdin_dict = {**common_dict}
    @property
    def name(self):
        return 'UNDEFINED _AMBERMover'
    @property
    def working_dir(self):
        '''str or None: Working directory in which intermediate files are
        stored. If set to None, will use a temporary directory that gets deleted
        once the AMBER operation is finished (not recommended). Current directory
        by default.'''
        return self._working_dir
    @working_dir.setter
    def working_dir(self, value):
        self._working_dir = value
    @property
    def prefix(self):
        return self._prefix
    @prefix.setter
    def prefix(self, value):
        self._prefix = value
    @property
    def solvent(self):
        '''Solvent or None: The type of solvent to use for an explicit-solvent
        simulation, if any. If this is set to None, an implicit-solvent
        simulation will be performed instead. None by default.'''
        return self._solvent
    @solvent.setter
    def solvent(self, value):
        self._solvent = value
        wet = value is not None
        self._min_mdin_dict['igb'] = int(not wet)
        self._min_mdin_dict['ntb'] = int(wet)
        self._mdin_dict['igb'] = int(not wet)
        self._mdin_dict['ntb'] = int(wet)
        self._amber_executable = 'pmemd' if wet else 'sander'
    @property
    def solvent_shape(self):
        '''SolventShape: The shape of solvent to use for an explicit-solvent
        simulation. ``OCT`` by default.'''
        return self._solvent_shape
    @solvent_shape.setter
    def solvent_shape(self, value):
        self._solvent_shape = value
    @property
    def cutoff(self):
        '''int or float: Van der Waals cutoff in angstroms. 15 by default.'''
        if self._min_mdin_dict['cut'] != self._mdin_dict['cut']:
            warnings.warn('Van der Waals cutoff values don\'t match for '
                          'minimization and simulation. To fix this, assign '
                          'something to cutoff. Returning simulation cutoff '
                          'value.',
                          RuntimeWarning)
        return self._mdin_dict['cut']
    @cutoff.setter
    def cutoff(self, value):
        self._min_mdin_dict['cut'] = value
        self._mdin_dict['cut'] = value
    @property
    def ref_pose(self):
        '''Pose or None: Reference Pose to give coordinate constraints, if any.
        It is safe to use the Pose that you are operating on for this purpose,
        if you want to constrain the Pose to its original coordinates; you just
        won't be able to use it again after an `apply`_, since the Pose will
        change.'''
        return self._ref_pose
    @ref_pose.setter
    def ref_pose(self, value):
        self._ref_pose = value
    @property
    def cst_mask(self):
        '''str or None: Atom mask string specifying which atoms are affected by
        the coordinate constraints. Given in ambmask syntax: see section 19.1 of
        the 2018 AMBER manual, on page 395.'''
        try:
            return self._mdin_dict['restraintmask']
        except KeyError:
            return None
    @cst_mask.setter
    def cst_mask(self, value):
        if value is None:
            try:
                del self._min_mdin_dict['restraintmask']
            except KeyError:
                pass
            try:
                del self._mdin_dict['restraintmask']
            except KeyError:
                pass
            return
        self._min_mdin_dict['restraintmask'] = value
        self._mdin_dict['restraintmask'] = value
    @property
    def cst_weight(self):
        '''int or float or None: Weight of coordinate constraints, in
        kcal/mol/angstroms^2.'''
        try:
            return self._mdin_dict['restraint_wt']
        except KeyError:
            return None
    @cst_weight.setter
    def cst_weight(self, value):
        if value is None:
            try:
                del self._min_mdin_dict['restraint_wt']
            except KeyError:
                pass
            try:
                del self._mdin_dict['restraint_wt']
            except KeyError:
                pass
            return
        self._min_mdin_dict['restraint_wt'] = value
        self._mdin_dict['restraint_wt'] = value
    @property
    def min_mdin_dict(self):
        '''Key-value pairs for AMBER minimization parameters. Only tamper with
        these if you've ever made a NAMELIST file for an AMBER simulation
        (because that's what this describes). For all available parameters, see
        section 17.6 of the 2018 AMBER manual, on page 320.'''

        return self._min_mdin_dict
    @min_mdin_dict.setter
    def min_mdin_dict(self, value):
        assert isinstance(value, dict), 'min_mdin_dict must be a dict.'
        self._min_mdin_dict = value
    @property
    def mdin_dict(self):
        '''dict: Key-value pairs for AMBER simulation parameters. Only tamper
        with these if you've ever made a NAMELIST file for an AMBER simulation
        (because that's what this describes). For all available parameters, see
        section 17.6 of the 2018 AMBER manual, on page 320.'''

        return self._mdin_dict
    @mdin_dict.setter
    def mdin_dict(self, value):
        assert isinstance(value, dict), 'mdin_dict must be a dict.'
        self._mdin_dict = value
    def _make_run_md_args(self, *, working_dir, prefix, minp):
        '''Create a dict of keyword arguments that can be passed to `run_md`_
        for a particular run, based on some local variables calculated during
        an `apply`_ call. Intended for internal use during said calls.

        Parameters
        ----------
        working_dir : str or TemporaryDirectory
            Object representing current working directory and behaving like a
            string.
        prefix : str
            Prefix used for this run.
        minp : bool
            Whether to produce args for a minimization (True) or a simulation
            (False).

        Returns
        -------
        dict
            Dict of arguments that can be unpacked into a `run_md`_ call during
            an `apply`_ call when all the necessary files already exist.
        '''

        suffix_dict = _AMBERMover.MIN_SUFFIXES if minp \
                      else _AMBERMover.MD_SUFFIXES
        args = {key: os.path.join(working_dir, prefix+suffix) \
                for key, suffix in suffix_dict.items()}
        args['executable'] = self._amber_executable
        args['overwrite'] = True
        if self.ref_pose is None:
            del args['refc']
        args['mdin_dict'] = copy.copy(self._min_mdin_dict if minp \
                                      else self._mdin_dict)
        return args
    def apply(self, pose):
        '''This half of apply will create a topology file for the pose, and
        minimize it in AMBER, which is the bare minimum for what you need to do
        before running a simulation. The topology file will be stored at
        ``[working_dir]/[prefix]-prmtop.parm7``, and the restart file for the
        minimized structure will be stored at
        ``[working_dir]/[prefix]-min-restrt.nc``. The local variables
        ``working_dir``, ``prefix``, and ``min_paths`` will be returned in a
        tuple, to be used for the second half, to be written in the subclass.

        Parameters
        ----------
        pose : rosetta.core.Pose
            Pose to minimize.

        Returns
        -------
        (str or TemporaryDirectory, str, dict)
            A tuple containing the local variables ``working_dir``,
            ``prefix``, and ``min_args``. ``working_dir`` is a string or
            TemporaryDirectory object representing the directory where the
            intermediate files are stored, which must be kept alive under most
            circumstances to prevent the directory from being deleted in the
            case that it's temporary. ``prefix`` is the string prefix for files
            created during this run; if the preperty `prefix`_ is set, it will
            be equal to it. ``min_args`` is a dict representing the arguments
            that were passed to `run_md`_ during the minimization that occurred
            when this `apply`_ was called.
        '''

        working_dir = self.working_dir if self.working_dir is not None \
                      else tempfile.TemporaryDirectory()
        prefix = _timestamp() if self.prefix is None else self.prefix
        leap_args = {key: os.path.join(working_dir, prefix+suffix) \
                     for key, suffix in _AMBERMover.LEAP_SUFFIXES.items()}
        leap_args['solvent'] = self.solvent
        leap_args['solvent_shape'] = self.solvent_shape
        pose_to_amber_params(pose, **leap_args)
        _check_pose_convertibility(pose,
                                   leap_args['crd_path'],
                                   leap_args['top_path'])
        min_args = self._make_run_md_args(working_dir=working_dir,
                                          prefix=prefix,
                                          minp=True)
        ref_pose = self.ref_pose
        if ref_pose is not None:
            ref_leap_args = copy.copy(leap_args)
            ref_leap_args['crd_path'] = min_args['refc']
            pose_to_amber_params(ref_pose, **leap_args)
            # Check whether AMBER-ified Pose has the same topology as the
            # original Pose. If it doesn't, then all our work would be for
            # naught (since leaving the missing residues in their places
            # wouldn't make sense, and deleting them would be rude), so we
            # raise a TopologyError and call it a day.
            _check_pose_convertibility(ref_pose,
                                       ref_leap_args['crd_path'],
                                       ref_leap_args['top_path'])
        run_md(**min_args)
        # Actual .apply methods aren't supposed to return anything, but this is
        # only half of one; all the local variables need to get carried over,
        # especially working_dir (so that it doesn't get destroyed in the case
        # that it's a temporary directory).
        return working_dir, prefix, min_args

## Classes supporting pose_to_amber_params

class Solvent(enum.Enum):
    '''A solvent choice for anything that requires you to set or provide a
    solvent, like `pose_to_amber_params`_ or `AMBERSimulateMover`_. Can be:

    SPCE_WATER
        SPC/E (extended simple point charge) waters.
    TIP3P_WATER
        TIP3P (transferable intramolecular potential with 3 points) waters.
    METHANOL
        Methanol.
    '''

    SPCE_WATER = 'SPCE_WATER'
    TIP3P_WATER = 'TIP3P_WATER'
    METHANOL = 'METHANOL'

class SolventShape(enum.Enum):
    '''The choice of shape for the periodic boundaries you're using for an
    explicit-solvent simulation, for anything that requires you to set or
    provide it, like `pose_to_amber_params`_ or `AMBERSimulateMover`_. Can be:

    BOX
        Rectangular prism.
    OCT
        Truncated octahedron.
    '''

    BOX = 'BOX'
    OCT = 'OCT'

### Front end classes

## Errors

class TopologyError(RuntimeError):
    '''An error encountered when there's something wrong with a topology.'''

class TopologySizeError(TopologyError):
    '''An error encountered when a topology has the wrong number of residues or
    atoms.'''

class TopologyTypeError(TopologyError):
    '''An error encountered when a topology has the wrong kind of residues.'''

## TrajToPoses

class TrajToPoses:
    '''A container that provides a ``Pose`` for each frame of a given pytraj
    ``Trajectory``. The Poses are created dynamically when you ask for them, so
    this object is not very memory-intensive. The trade-off, however, is that
    getting a ``Pose`` from it is fairly expensive, so try not to do it more
    times than you need to. (One could imagine a version of this object that
    stores weak references to all the ``Pose``s it emits, so that it returns the
    same ``Pose`` object every time for a given index unless it stops being used
    elsewhere. But then you'd have to make a new TrajToPoses_ object
    every time you wanted the original ``Pose`` corresponding to that index,
    which is unidiomatic.)

    Example
    -------
    This is how one would obtain a ``Pose`` for the first frame of a
    ``Trajectory``::

        import pyrosetta as pr; pr.init()
        import pytraj as pt
        from miniambrose import TrajToPose

        traj = pt.iterload('my-traj.mdcrd', 'my-top.prmtop')
        poses = TrajToPose(traj)
        pose = poses[0]

    Parameters
    ----------
    traj : pytraj.Trajectory
        The trajectory the Poses are derived from.
    options : rosetta.core.import_pose.ImportPoseOptions, optional
        The options used to derive the Poses.
    __maps : tuple(_ResidueMap), optional
        A tuple of `_ResidueMap`_s to use as its internal map of Rosetta
        residues to pytraj residues. Used whenever slices of this object are
        made.

    Attributes
    ----------
    traj : pytraj.Trajectory
        The trajectory the Poses are derived from.
    options : rosetta.core.import_pose.ImportPoseOptions
        The options used to derive the Poses.
    __maps : tuple(_ResidueMap)
        A tuple of `_ResidueMap`_s to use as its internal map of Rosetta
        residues to pytraj residues.
    '''

    def __init__(self, traj, options=None, *, __maps=None):
        self.traj = traj
        self.options = options or pr.rosetta.core.import_pose.ImportPoseOptions
        if __maps is not None:
            self.__maps = __maps
        else:
            self.__maps = _TopologyParser().parse_topology(traj.top)
    def __len__(self):
        return len(self.traj)
    def __iter__(self):
        for i in range(len(self)):
            yield self[i]
    def __getitem__(self, key):
        if isinstance(key, slice):
            return TrajToPoses(self.traj[key], __maps=self.__maps)

        frame = self.traj[key]
        sfr  = pr.rosetta.core.io.StructFileRep()
        ChainAtoms      = pr.rosetta.utility.vector0_core_io_AtomInformation
        chain_atoms_map = collections.defaultdict(ChainAtoms)
        atom_i = itertools.count(1)
        for res_i, res_map in enumerate(self.__maps, 1):
            for atom_record in res_map:
                if atom_record.pt_i is None:
                    continue
                ai = atom_record.atom_information
                ai.serial = next(atom_i)
                ai.resSeq = res_i
                ai.x, ai.y, ai.z = frame[atom_record.pt_i]
                chain_atoms_map[ai.chainID].append(ai)
        for chain in chain_atoms_map.values():
            sfr.chains().append(chain)
        sfr.header().finalize_parse()

        pose = pr.Pose()
        rts = _global_residue_type_set()
        options = pr.rosetta.core.import_pose.ImportPoseOptions()
        pr.rosetta.core.import_pose.build_pose_as_is(sfr, pose, rts, options)
        return pose

class AMBERMinMover(_AMBERMover):
    '''A Mover-like object that uses sander or PMEMD to minimize a pose.

    Intermediate files are stored by default in the working directory, prefixed
    with a timestamp. The directory and prefix may be overridden. Output of the
    files can be turned off altogether by setting the directory to ``None``,
    although the files will still be output intermediately and cleaned up at the
    end. Here is a list of all files output by this mover:

    <prefix>-leap.log
        The log file where the LEaP output for the conversion of the pose(s)
        to AMBER input files is recorded.
    <prefix>-prmtop.parm7
        The topology file describing the chemistry and connectivity of the
        protein(s) in your minimization, generated by LEaP.
    <prefix>-inpcrd.rst7
        The initial coordinates for your minimization, generated by LEaP.
    <prefix>-min-mdin.in
        The FORTRAN NAMELIST file that contains the minimization parameters.
    <prefix>-min-mdout.out
        The log file where each ``mdinfo`` is recorded during the minimization.
    <prefix>-min-mdinfo
        The log file where high-level information about the current state of the
        minimization is recorded.
    <prefix>-min-restrt.nc
        The latest coordinates for your minimization. After the minimization is
        over, this gives the new coordinates for your input Pose.

    The following file may additionally be output if a reference Pose (with
    necessarily the same topology as your simulated Pose) is specified:

    <prefix>-refc.rst7
        The coordinates of the reference structure, generated by LEaP.

    Member variables are Python properties, as opposed to being accessed via
    getters and setters as in a traditional Mover. See the getters for their
    documentation, or see the constructor parameters given below, which
    correspond to member variables.

    This object provides, but does not implement or document, all the methods of
    the base Mover class. The only method from the base class that actually
    works is `apply`_.

    Parameters
    ----------
    working_dir : str or None, optional
        The directory to place the intermediate files in. Defaults to the
        current working directory. If set to None, intermediate files will
        be output to a temporary directory, and deleted once the mover has
        finished running. (This is not recommended, especially if you're running
        the mover multiple times on the same structure.)
    prefix : str or None, optional
        The prefix of the intermediate files output by the simulation. If set to
        None, will use the timestamp at which `apply`_ was called.
    solvent : Solvent or None
        The type of solvent to use for an explicit-solvent simulation, if any.
        If this is set to None, an implicit-solvent simulation will be performed
        instead. None by default.
    solvent_shape : SolventShape
        The shape of solvent to use for an explicit-solvent simulation. ``OCT``
        by default.
    cutoff : int or float, optional
        Van der Waals cutoff in angstroms. 15 by default.
    ref_pose : rosetta.core.Pose, optional
        Reference pose to give coordinate constraints, if any. It is safe to
        use the pose that you are operating on for this purpose, if you want
        to constrain the pose to its original coordinates.
    cst_mask : str, optional
        Atom mask specifying which atoms are affected by the coordinate
        constraints. Given in ambmask syntax; see section 19.1 of the 2018
        AMBER manual, on page 395.
    cst_weight : float, optional
        Weight of coordinate constraints, in kcal/mol/angstroms^2.
    '''

    def __init__(self, working_dir='', prefix=None, solvent=None,
                 solvent_shape=SolventShape.OCT, cutoff=15.):
        _AMBERMover.__init__(self,
                             working_dir=working_dir,
                             prefix=prefix,
                             solvent=solvent,
                             solvent_shape=solvent_shape,
                             cutoff=cutoff)
    @property
    def name(self):
        return 'AMBERMinMover'
    def apply(self, pose):
        '''Minimize the given pose using sander or PMEMD.

        Parameters
        ----------
        pose : rosetta.core.Pose
            Pose to minimize.
        '''

        # pylint: disable=unused-variable
        # (working_dir needs to be alive so that the directory can continue
        #  existing if it's a temporary directory)
        working_dir, _, min_args = _AMBERMover.apply(self, pose)
        _transfer_xyz_from_pose_to_pose(
            TrajToPoses(pt.iterload(min_args['restrt'],
                                    min_args['prmtop']))[0],
            pose)
        # working_dir should get cleaned up on its own if it's a
        # TemporaryDirectory

class AMBERSimulateMover(_AMBERMover):
    '''A Mover-like object that uses sander or PMEMD to minimize and then
    simulate a pose, ultimately overwriting it with a frame from the
    simulation's trajectory.

    Intermediate files are stored by default in the working directory, prefixed
    with a timestamp. The directory and prefix may be overridden. Output of the
    files can be turned off altogether by setting the directory to ``None``,
    although the files will still be output intermediately and cleaned up at the
    end. Here is a list of all files output by this mover:

    <prefix>-leap.log
        The log file where the LEaP output for the conversion of the pose(s)
        to AMBER input files is recorded.
    <prefix>-prmtop.parm7
        The topology file describing the chemistry and connectivity of the
        protein(s) in your simulation, generated by LEaP.
    <prefix>-inpcrd.rst7
        The initial coordinates for your minimization, generated by LEaP.
    <prefix>-min-mdin.in
        The FORTRAN NAMELIST file that contains the minimization parameters.
    <prefix>-min-mdout.out
        The log file where each ``mdinfo`` is recorded during the minimization.
    <prefix>-min-mdinfo
        The log file where high-level information about the current state of the
        minimization is recorded.
    <prefix>-min-restrt.nc
        The latest coordinates for your minimization. After the minimization is
        over, this gives the starting coordinates for your simulation.
    <prefix>-mdin.in
        The FORTRAN NAMELIST file that contains the simulation parameters.
    <prefix>-mdout.out
        The log file where each ``mdinfo`` is recorded during the simulation.
    <prefix>-mdinfo
        The log file where high-level information about the current state of the
        simulation is recorded.
    <prefix>-restrt.nc
        The latest coordinates and velocities for your simulation.
    <prefix>-mdcrd.nc
        The frames for your simulation's trajectory. Your Pose's new
        conformation will be selected from among these by `pose_selector`_.

    The following file may additionally be output if a reference Pose (with
    necessarily the same topology as your simulated Pose) is specified:

    <prefix>-refc.rst7
        The coordinates of the reference structure, generated by LEaP.

    Member variables are Python properties, as opposed to being accessed via
    getters and setters as in a traditional Mover. See the getters for their
    documentation, or see the constructor parameters given below, which
    correspond to member variables. The only writable member variable that does
    not appear in the constructor as a parameter is `pose_selector`_, which is
    the function that determines which frame to use from the simulated
    trajectory to overwrite the input Pose.

    This object provides, but does not implement or document, all the methods of
    the base Mover class. The only method from the base class that actually
    works is `apply`_.

    Parameters
    ----------
    working_dir : str or None, optional
        The directory to place the intermediate files in. Defaults to the
        current working directory. If set to None, intermediate files will
        be output to a temporary directory, and deleted once the mover has
        finished running. (This is not recommended, especially if you're running
        the mover multiple times on the same structure.)
    prefix : str or None, optional
        The prefix of the intermediate files output by the simulation. If set to
        None, will use the timestamp at which `apply`_ was called.
    solvent : Solvent or None
        The type of solvent to use for an explicit-solvent simulation, if any.
        If this is set to None, an implicit-solvent simulation will be performed
        instead. None by default.
    solvent_shape : SolventShape
        The shape of solvent to use for an explicit-solvent simulation. ``OCT``
        by default.
    duration : int or float, optional
        Duration of simulation in picoseconds. Must be set before the simulation
        runs.
    temperature : int or float, optional
        Temperature of simulation in kelvins. Must be set before the simulation
        runs.
    starting_temperature : int or float, optional
        Starting temperature of simulation in kelvins. 0 by default.
    cutoff : int or float, optional
        Van der Waals cutoff in angstroms. 15 by default.
    output_interval : int or float, optional
        Amount of time to wait in between writing frames and log data, in
        picoseconds. 1 by default.
    seed : int, optional
        Seed to use for random operations (e.g. Langevin thermostat).
        Nonnegative integers are valid seeds, while negative integers cause
        the seed to be determined by the hardware.
    ref_pose : rosetta.core.Pose, optional
        Reference pose to give coordinate constraints, if any. It is safe to
        use the pose that you are operating on for this purpose, if you want
        to constrain the pose to its original coordinates.
    cst_mask : str, optional
        Atom mask specifying which atoms are affected by the coordinate
        constraints. Given in ambmask syntax; see section 19.1 of the 2018
        AMBER manual, on page 395.
    cst_weight : float, optional
        Weight of coordinate constraints, in kcal/mol/angstroms^2.
    '''
    def __init__(self, working_dir='', prefix=None, solvent=None,
                 solvent_shape=SolventShape.OCT, duration=None,
                 temperature=None, starting_temperature=0, cutoff=15.,
                 output_interval=1., seed=-1, ref_pose=None, cst_mask=None,
                 cst_weight=None):
        _AMBERMover.__init__(self,
                             working_dir=working_dir,
                             prefix=prefix,
                             solvent=None,
                             solvent_shape=SolventShape.OCT,
                             cutoff=cutoff,
                             ref_pose=ref_pose,
                             cst_mask=cst_mask,
                             cst_weight=cst_weight)
        self._duration_is_int = isinstance(duration, int)
        if duration is not None:
            self.duration = duration
        if temperature is not None:
            self.temperature = temperature
        self.starting_temperature = starting_temperature
        self._output_interval_is_int = isinstance(output_interval, int)
        self.output_interval = output_interval
        self.seed = seed
        self._pose_selector = pose_selectors.last
    @property
    def name(self):
        return 'AMBERSimulateMover'
    @property
    def pose_selector(self):
        '''A function that selects which frame of the simulated trajectory to
        use to overwrite your input Pose when `apply`_ is called. In practice,
        it can be any function that takes an ordered list of Poses and returns
        one of them, since it operates on a TrajToPose_ object and returns one
        of its entries. By default, it's equal to the function `last`_ from the
        submodule `pose_selectors`_, which selects the last frame of the
        trajectory; there is also the more sophisticated option of
        `lowest_energy_from_end`_ available from the same module.'''
        return self._pose_selector
    @pose_selector.setter
    def pose_selector(self, value):
        self._pose_selector = value
    @property
    def duration(self):
        '''int or float: Duration of simulation in picoseconds.'''
        if self._duration_is_int:
            return self._mdin_dict['nstlim'] // 1000
        return self._mdin_dict['nstlim'] / 1000.
    @duration.setter
    def duration(self, value):
        self._duration_is_int = isinstance(value, int)
        self._mdin_dict['nstlim'] = int(value * 1000)
    @property
    def temperature(self):
        '''int or float: Temperature of simulation in kelvins.'''
        return self._mdin_dict['temp0']
    @temperature.setter
    def temperature(self, value):
        self._mdin_dict['temp0'] = value
    @property
    def starting_temperature(self):
        '''int or float: Starting temperature of simulation in kelvins. 0 by
        default.'''
        return self._mdin_dict['tempi']
    @starting_temperature.setter
    def starting_temperature(self, value):
        self._mdin_dict['tempi'] = value
    @property
    def output_interval(self):
        '''int or float: Amount of in-simulation time to wait in between writing
        successive frames and log data, in picoseconds. 1 by default.'''
        if self._output_interval_is_int:
            return self._mdin_dict['ntwx'] // 1000
        return self._mdin_dict['ntwx'] / 1000.
    @output_interval.setter
    def output_interval(self, value):
        self._output_interval_is_int = isinstance(value, int)
        self._mdin_dict['ntwx'] = int(value * 1000)
        self._mdin_dict['ntpr'] = int(value * 1000)
    @property
    def seed(self):
        '''int: Seed to use for random operations (in our case, the Langevin
        thermostat, if it's a wet simulation). Nonnegative integers are valid
        seeds, while negative integers cause the seed to be determined by the
        hardware. -1 by default.'''
        try:
            return self._mdin_dict['ig']
        except KeyError:
            return None
    @seed.setter
    def seed(self, value):
        if value is None:
            try:
                del self._mdin_dict['ig']
            except KeyError:
                return
        self._mdin_dict['ig'] = value
    def apply(self, pose):
        '''Minimize, then simulate the given pose using sander or PMEMD.

        Parameters
        ----------
        pose : rosetta.core.Pose
            Pose to minimize and simulate.
        '''

        assert self.duration is not None, 'Set duration first.'
        assert self.temperature is not None, 'Set temperature first.'

        # pylint: disable=unused-variable
        # (working_dir needs to be alive so that the directory can continue
        #  existing if it's a temporary directory)
        working_dir, prefix, _ = _AMBERMover.apply(self, pose)
        md_args = self._make_run_md_args(working_dir=working_dir,
                                         prefix=prefix,
                                         minp=False)
        run_md(**md_args)
        _transfer_xyz_from_pose_to_pose(
            self.pose_selector( # pylint: disable=not-callable
                TrajToPoses(
                    pt.iterload(md_args['mdcrd'],
                                md_args['prmtop']))),
            pose)
        # working_dir should get cleaned up on its own if it's a
        # TemporaryDirectory

def pose_to_amber_params(pose, crd_path, top_path, log_path='leap.log', *,
                         leap_script_dump_path=None, solvent=None,
                         solvent_shape=SolventShape.OCT, add_ions=True):
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
    log_path : str
        Path to the output LEaP log file. Default is 'leap.log' in the current
        directory.
    leap_script_dump_path : str or None
        Path to dump the generated LEaP script at. Default is None.
    solvent : Solvent or None
        What explicit solvent to use, if any. See the `Solvent`_ class for your
        options. None by default.
    solvent_shape : SolventShape
        What shape to use for the explicit solvent box, if any. See the
        `SolventShape`_ class for your options. ``OCT`` by default.
    add_ions : bool
        Whether to add neutralizing sodium or chloride ions in the case that
        your simulation is explicit-solvent and the charge of your protein is
        nonzero. True by default.
    '''

    if leap_script_dump_path is not None:
        leap_in_fd, leap_in_path = leap_script_dump_path, leap_script_dump_path
    else:
        leap_in_fd, leap_in_path = tempfile.mkstemp()
    try:
        with tempfile.NamedTemporaryFile() as pdb_f:
            ## dump the PDB we're going to feed to LEaP
            dump_amber_pdb(pose, pdb_f.name)

            ## generate and write the LEaP instructions
            leap_command_str = bytearray(
                open(os.path.join(AMBROSE_DIR, 'pose-to-traj.in.prototype')).read(),
                'ascii')
            # To ensure that this is not a bottleneck, we give the indices of the
            # replacement points in advance. We could technically forego the names
            # to gain even more time, but then it becomes pretty unreadable.
            if solvent:
                solvent_box, solvent_cmd = \
                    pose_to_amber_params.SOLVENTS[solvent]
            replacements = {('+LOG-FILE+', (8,)):
                                log_path,
                            ('+LOAD-SOLVENT+', (48,)):
                                solvent_cmd if solvent else '',
                            ('+PDB-PATH+', (81,)):
                                pdb_f.name,
                            ('+SOLVATEP+', (92,)):
                                '' if solvent else '#',
                            ('+SOLVENT-SHAPE+', (110,)):
                                solvent_shape.lower(),
                            ('+SOLVENT+', (133,)):
                                solvent_box if solvent else '',
                            ('+ADD-IONS-P+', (148, 182)):
                                '' if add_ions and solvent else '#',
                            ('+OUT-TOP-PATH+', (237,)):
                                top_path,
                            ('+OUT-CRD-PATH+', (252,)):
                                crd_path}
            # cumulative difference between given index and actual replacement
            # index, from accumulated differences in lengths between names and
            # their replacements:
            cum_diff = 0
            for (name, indices), replacement in replacements.items():
                len_name = len(bytearray(name, 'ascii'))
                byte_replacement = bytearray(replacement, 'ascii')
                for index in indices:
                    leap_command_str[index+cum_diff:index+len_name+cum_diff] = \
                        byte_replacement
                    cum_diff += len(byte_replacement) - len_name
            open(leap_in_fd, 'wb').write(leap_command_str)

            ## get LEaP to write the file
            subprocess.run([os.path.join(amber_bin(), 'tleap'),
                            '-f', leap_in_path])
    finally:
        if leap_script_dump_path is None:
            os.remove(leap_in_path)
## the "load" commands for different solvents, used to construct the LEaP
## instructions
pose_to_amber_params.SOLVENTS = \
    {Solvent.SPCE_WATER:
     ('SPCBOX',
      'source leaprc.water.spce'),
     Solvent.TIP3P_WATER:
     ('TIP3PBOX',
      'source leaprc.water.tip3p'),
     Solvent.METHANOL:
     ('MEOHBOX',
      'loadoff solvents.lib\nloadamberparams frcmod.meoh')}
pose_to_amber_params.SOLVENT_SHAPES = \
    {Solvent.BOX: 'box',
     Solvent.OCT: 'oct'}

def run_md(executable='sander', *, overwrite=False, mdin=None,
           mdin_dict=None, mdout=None, prmtop=None, inpcrd=None, restrt=None,
           refc=None, mdcrd=None, mdinfo=None):
    '''Runs an AMBER MD program (sander by default) with the given arguments, with
    the additional feature that the mdin file, which normally gives the
    simulation parameters, can be replaced with a dict specifying the same
    parameters. All arguments are optional, as in the sander executable; if a
    value is needed but not specified, sander will assume the default (which is
    the argument's name in this case).

    If you don't want to make the mdin_dict parameters by hand, the templates_
    submodule can synthesize some common formulas for you.

    Parameters
    ----------
    executable : str, optional
        The AMBER executable to run the MD with. sander by default.
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
        result = subprocess.run([os.path.join(amber_bin(), executable)] + arguments)
        _debug_print(f'{executable} run with result: {result}')
        return result
