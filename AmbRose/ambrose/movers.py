'''The _AMBERMovers and their supporting functions and classes. These fake
movers all convert a Pose to AMBER params, do something with them in AMBER
(either a minimization or simulation, so far), and replace the Pose with
the result in some way.'''

import datetime
import weakref
import copy
import tempfile
import os
import sys
import warnings

# pragma pylint: disable=import-error
import pytraj as pt
# pragma pylint: enable=import-error

from . import utils
from . import enums
from . import validation
from . import traj_to_poses
from . import pose_to_traj
from . import pose_selectors

### Supporting functions

def _timestamp():
    '''Gets timestamp in standardized format.'''

    return datetime.datetime.now().isoformat()

def _transfer_xyz_from_pose_to_pose(source_pose, dest_pose):
    '''Copy coordinates from one pose to another, assuming that they have the
    exact same topology (we don't bother to check).'''

    for i in range(1, dest_pose.size()+1):
        dest_r = dest_pose.residue(i)
        source_r = source_pose.residue(i) # pylint: disable=no-member
        for j in range(1, dest_r.natoms()):
            dest_r.set_xyz(j, source_r.xyz(j))

### Classes

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
                 solvent_shape=enums.SolventShape.OCT, cutoff=15.,
                 ref_pose=None, cst_mask=None, cst_weight=None):
        self._working_dir = working_dir
        self._prefix = prefix
        self._ref_pose = ref_pose
        self._solvent = solvent
        self._solvent_shape = solvent_shape
        wet = solvent is not None
        self._amber_executable = 'pmemd.cuda'
        common_dict = {'igb': int(not wet)*8,
                       'ntb': int(wet),
                       'cut': cutoff}
        if ref_pose is not None:
            common_dict['ntr'] = 1
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
        '''str or None: The prefix of the intermediate files output by the
        simulation. If set to None, will use the timestamp at which `apply`_
        was called.'''
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
        self._min_mdin_dict['igb'] = int(not wet)*8
        self._min_mdin_dict['ntb'] = int(wet)
        self._mdin_dict['igb'] = int(not wet)*8
        self._mdin_dict['ntb'] = int(wet)
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
        '''int or float: Van der Waals cutoff in angstroms. Ignored for
        implicit-solvent simulations, since pmemd.cuda doesn't allow a cutoff
        for those. 15 by default.'''
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
        if value is None:
            if 'ntr' in self._min_mdin_dict:
                del self._min_mdin_dict['ntr']
            if 'ntr' in self._mdin_dict:
                del self._mdin_dict['ntr']
        else:
            self._min_mdin_dict['ntr'] = 1
            self._mdin_dict['ntr'] = 1
    @property
    def cst_mask(self):
        '''str or None: Atom mask string specifying which atoms are affected by
        the coordinate constraints. Given in ambmask syntax: see section 19.1 of
        the 2019 AMBER manual, on page 410.'''
        return self._mdin_dict.get('restraintmask')
    @cst_mask.setter
    def cst_mask(self, value):
        if value is None:
            if 'restraintmask' in self._min_mdin_dict:
                del self._min_mdin_dict['restraintmask']
            if 'restraintmask' in self._mdin_dict:
                del self._mdin_dict['restraintmask']
        else:
            self._min_mdin_dict['restraintmask'] = value
            self._mdin_dict['restraintmask'] = value
    @property
    def cst_weight(self):
        '''int or float or None: Weight of coordinate constraints, in
        kcal/mol/angstroms^2.'''
        return self._mdin_dict.get('restraint_wt')
    @cst_weight.setter
    def cst_weight(self, value):
        if value is None:
            if 'restraint_wt' in self._min_mdin_dict:
                del self._min_mdin_dict['restraint_wt']
            if 'restraint_wt' in self._mdin_dict:
                del self._mdin_dict['restraint_wt']
        else:
            self._min_mdin_dict['restraint_wt'] = value
            self._mdin_dict['restraint_wt'] = value
    @property
    def min_mdin_dict(self):
        '''dict: Key-value pairs for AMBER minimization parameters. Only tamper
        with these if you've ever made a NAMELIST file for an AMBER simulation
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
        # Necessary so that pmemd.cuda doesn't barf. Why don't we just leave
        # out the cutoff altogether? Because, apparently, the default of 0 is
        # not good enough, because it's less than 999 angstroms. You'd think
        # that this means that the cutoff must be large but finite, but no,
        # anything above 999 angstroms gets treated as infinite anyway. Aaaargh.
        if args['mdin_dict']['igb'] != 0:
            args['mdin_dict']['cut'] = 1000.
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

        if self.working_dir is not None:
            assert os.path.exists(self.working_dir) or not self.working_dir, \
                   'Working directory must exist.'
        working_dir = self.working_dir if self.working_dir is not None \
                      else tempfile.TemporaryDirectory()
        prefix = _timestamp() if self.prefix is None else self.prefix
        leap_args = {key: os.path.join(working_dir, prefix+suffix) \
                     for key, suffix in _AMBERMover.LEAP_SUFFIXES.items()}
        leap_args['solvent'] = self.solvent
        leap_args['solvent_shape'] = self.solvent_shape
        pose_to_traj.pose_to_amber_params(pose, **leap_args)
        validation.check_file_existence(leap_args['crd_path'],
                                        'input coordinates file')
        validation.check_file_existence(leap_args['top_path'],
                                        'input topology file')
        validation.check_pose_convertibility(
            pose,
            leap_args['crd_path'],
            leap_args['top_path'])
        min_args = self._make_run_md_args(working_dir=working_dir,
                                          prefix=prefix,
                                          minp=True)
        ref_pose = self.ref_pose
        if ref_pose is not None:
            ref_leap_args = copy.copy(leap_args)
            ref_leap_args['crd_path'] = min_args['refc']
            pose_to_traj.pose_to_amber_params(ref_pose, **ref_leap_args)
            validation.check_pose_convertibility(
                ref_pose,
                ref_leap_args['crd_path'],
                ref_leap_args['top_path'])
        pose_to_traj.run_md(**min_args)
        # Actual .apply methods aren't supposed to return anything, but this is
        # only half of one; all the local variables need to get carried over,
        # especially working_dir (so that it doesn't get destroyed in the case
        # that it's a temporary directory).
        return working_dir, prefix, min_args

### Premier classes

class AMBERMinMover(_AMBERMover):
    '''A Mover-like object that uses pmemd.cuda to minimize a pose.

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
    solvent : Solvent or None, optional
        The type of solvent to use for an explicit-solvent simulation, if any.
        If this is set to None, an implicit-solvent simulation will be performed
        instead. None by default.
    solvent_shape : SolventShape, optional
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
        constraints. Given in ambmask syntax; see section 19.1 of the 2019
        AMBER manual, on page 410.
    cst_weight : float, optional
        Weight of coordinate constraints, in kcal/mol/angstroms^2.
    '''

    def __init__(self, working_dir='', prefix=None, solvent=None,
                 solvent_shape=enums.SolventShape.OCT, cutoff=15.):
        _AMBERMover.__init__(self,
                             working_dir=working_dir,
                             prefix=prefix,
                             solvent=solvent,
                             solvent_shape=solvent_shape,
                             cutoff=cutoff)
    def __repr__(self):
        return f'AMBERMinMover(\n' \
               + f'  working_dir={self.working_dir},\n' \
               + f'  prefix={self.prefix},\n' \
               + f'  solvent={self.solvent},\n' \
               + f'  solvent_shape={self.solvent_shape},\n' \
               + f'  cutoff={self.cutoff})'
    @property
    def name(self):
        return 'AMBERMinMover'
    def apply(self, pose):
        '''Minimize the given pose using pmemd.cuda.

        Parameters
        ----------
        pose : rosetta.core.Pose
            Pose to minimize.
        '''

        # pylint: disable=unused-variable
        # (working_dir needs to be alive so that the directory can continue
        #  existing if it's a temporary directory)
        working_dir, _, min_args = _AMBERMover.apply(self, pose)
        validation.check_file_existence(
            min_args['restrt'],
            'minimization results')
        _transfer_xyz_from_pose_to_pose(
            traj_to_poses.pose_from_amber_params(
                min_args['restrt'],
                min_args['prmtop']),
            pose)
        # working_dir should get cleaned up on its own if it's a
        # TemporaryDirectory

class AMBERSimulateMover(_AMBERMover):
    '''A Mover-like object that uses pmemd.cuda to minimize and then simulate a
    pose, ultimately overwriting it with a frame from the simulation's
    trajectory.

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

    This object provides, but does not implement or document, all the methods
    of the base Mover class. See the documentation for the class `_NotAMover`_
    for more info.

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
    solvent : Solvent or None, optional
        The type of solvent to use for an explicit-solvent simulation, if any.
        If this is set to None, an implicit-solvent simulation will be performed
        instead. None by default.
    solvent_shape : SolventShape, optional
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
        constraints. Given in ambmask syntax; see section 19.1 of the 2019
        AMBER manual, on page 410.
    cst_weight : float, optional
        Weight of coordinate constraints, in kcal/mol/angstroms^2.
    '''

    def __init__(self, working_dir='', prefix=None, solvent=None,
                 solvent_shape=enums.SolventShape.OCT, duration=None,
                 temperature=None, starting_temperature=0, cutoff=15.,
                 output_interval=1., seed=-1, ref_pose=None, cst_mask=None,
                 cst_weight=None):
        _AMBERMover.__init__(
            self,
            working_dir=working_dir,
            prefix=prefix,
            solvent=solvent,
            solvent_shape=solvent_shape,
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
    def __repr__(self):
        return f'AMBERMinMover(\n' \
               + f'  working_dir={self.working_dir},\n' \
               + f'  prefix={self.prefix},\n' \
               + f'  solvent={self.solvent},\n' \
               + f'  solvent_shape={self.solvent_shape},\n' \
               + f'  duration={self.duration},\n' \
               + f'  temperature={self.temperature},\n' \
               + f'  starting_temperature={self.starting_temperature},\n' \
               + f'  cutoff={self.cutoff},\n' \
               + f'  output_interval={self.output_interval},\n' \
               + f'  seed={self.seed},\n' \
               + f'  pose_selector={self.pose_selector},\n' \
               + f'  ref_pose={self.ref_pose},\n' \
               + f'  cst_mask={self.cst_mask},\n' \
               + f'  cst_weight={self.cst_weight})'
    @property
    def name(self):
        return 'AMBERSimulateMover'
    @property
    def pose_selector(self):
        '''function: A function that selects which frame of the simulated
        trajectory to use to overwrite your input Pose when `apply`_ is called.
        In practice, it can be any function that takes an ordered list of Poses
        and returns one of them, since it operates on a `TrajToPoses`_ object
        and returns one of its entries. By default, it's equal to the function
        `last`_ from the submodule `pose_selectors`_, which selects the last
        frame of the trajectory; there is also the more sophisticated option of
        `lowest_energy_from_end`_ available from the same module.'''
        return self._pose_selector
    @pose_selector.setter
    def pose_selector(self, value):
        self._pose_selector = value
    @property
    def duration(self):
        '''int or float: Duration of simulation in picoseconds.'''
        if 'nstlim' not in self._mdin_dict:
            return None
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
        if 'temp0' not in self._mdin_dict:
            return None
        return self._mdin_dict.get('temp0')
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
    def simulate(self, pose):
        '''Minimize, then simulate the given Pose using pmemd.cuda. Return a
        Trajectory object of the frames output during the simulation.

        Parameters
        ----------
        pose : rosetta.core.Pose
            Pose to minimize and simulate.

        Returns
        -------
        pytraj.Trajectory
            The frames from the simulation.
        '''

        assert self.duration is not None, 'Set duration first.'
        assert self.temperature is not None, 'Set temperature first.'

        # pylint: disable=unused-variable
        # (working_dir needs to be alive so that the directory can continue
        #  existing if it's a temporary directory)
        working_dir, prefix, min_args = _AMBERMover.apply(self, pose)
        validation.check_file_existence(
            min_args['restrt'],
            'minimization results')
        md_args = self._make_run_md_args(working_dir=working_dir,
                                         prefix=prefix,
                                         minp=False)
        pose_to_traj.run_md(**md_args)
        return traj_to_poses.TrajToPoses(pt.iterload(md_args['mdcrd'],
                                                     md_args['prmtop']))
        # working_dir should get cleaned up on its own if it's a
        # TemporaryDirectory
    def apply(self, pose):
        '''Minimize, then simulate the given pose using pmemd.cuda. Replace the
        Pose with a frame from the simulation selected by this mover's
        `pose_selector`_.

        Parameters
        ----------
        pose : rosetta.core.Pose
            Pose to minimize and simulate.
        '''

        _transfer_xyz_from_pose_to_pose(
            self.pose_selector( # pylint: disable=not-callable
                self.simulate(pose)),
            pose)
