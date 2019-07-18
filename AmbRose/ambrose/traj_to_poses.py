'''The `TrajToPoses`_ class and its supporting functions.'''

# pragma pylint: disable=bad-whitespace

import collections
import itertools

from functools import reduce

# pragma pylint: disable=import-error
import pyrosetta as pr
import pytraj as pt
# pragma pylint: enable=import-error

from . import utils

### Supporting functions

def _chemical_manager():
    return pr.rosetta.core.chemical.ChemicalManager.get_instance()

# Currently an orphan, but may be needed for future expansions of
# _TopologyParser.
def _global_element_set():
    if not hasattr(_global_element_set, 'es'):
        _global_element_set.es = _chemical_manager().element_set('default')
    return _global_element_set.es

def _global_residue_type_set():
    if not hasattr(_global_residue_type_set, 'rts'):
        _global_residue_type_set.rts = \
            _chemical_manager().residue_type_set('fa_standard')
    return _global_residue_type_set.rts

def _get_pr_res(name):
    return _global_residue_type_set().name_map(name)

def _patch_pr_res(res, patch_name):
    return _global_residue_type_set().patch_map()[patch_name][1].apply(res)

### Supporting classes

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
        ai = pr.rosetta.core.io.AtomInformation()
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
    chain_id : str, optional
        A single printable ASCII character representing the chain to which this
        residue belongs. '@' by default (0x40).

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
        residue belongs.
    name3 : str
        The three-letter residue name for this residue in Rosetta.
    atoms : list(:obj:`_AtomRecord`)
        A list of the :obj:`_AtomRecord`s of this residue.
    __names_to_indices : dict
        A dict that maps the 4-character names of the atoms in the corresponding
        Rosetta ``ResidueType`` to their indices in `atoms`.
    '''

    def __init__(self, name, patches, first_atom_i, always_base_p=False,
                 chain_id='@'):
        self.name = name
        self.patches = patches
        self.first_atom_i = first_atom_i
        self.chain_id = chain_id

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
        'CYM':'CYZ',
        'CYX':'CYD',
        'GLH':'GLU',
        'HID':'HIS_D',
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
             END:    {_TopologyParser.NAA: (PRE_AA, ('CHAIN_INCR', 'NAA')),
                      _TopologyParser.AA:  (PRE_AA, ('CHAIN_INCR', 'NAA')),
                      frozenset({'ACE'}):  (ACE,    ('CHAIN_INCR',)),
                      _TopologyParser.NA:  (NA,     ('CHAIN_INCR', 'NA')),
                      _TopologyParser.NA3: (NA,     ('CHAIN_INCR', 'NA3')),
                      _TopologyParser.NA5: (NA,     ('CHAIN_INCR', 'NA5'))}}
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
        chain_id = 'A'

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
            utils.debug_print('residue consumed:', pt_res_i, pt_res)

        for residue in top.residues:
            utils.debug_print('residue:', residue)
            tokens = list(self.transit(residue.name))
            if pt_res_i == top.n_residues-2 and tokens[-1] == 'AA':
                tokens.append('CAA')
            for token in tokens:
                utils.debug_print('token:', token)

                ## deal with "control" tokens:
                if token == 'SKIP':
                    consume_residue()
                    continue
                if token == 'CHAIN_INCR':
                    if ord(chain_id) < 126:
                        chain_id = chr(ord(chain_id) + 1)
                    continue

                naa_or_caa_p = token in ('NAA', 'CAA')
                na3_or_na5_p = token in ('NA3', 'NA5')

                ## determine name of new residue, name
                pt_name = pt_res.name
                name = _TopologyParser.PT_TO_PR_RES_NAMES.get(pt_name, pt_name)
                if naa_or_caa_p and len(name) == 4 and name[0] == 'N':
                    name = name[1:]

                ## determine patches to be applied, patcheshelp
                patches = _TopologyParser.TOKENS_TO_PATCHES.get(token, ())

                ## construct _ResidueMap records:
                res_map = _ResidueMap(name, patches, pt_atom_i,
                                      always_base_p=(naa_or_caa_p or \
                                                     na3_or_na5_p),
                                      chain_id=chain_id)
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

### The money class

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

### Front-facing functions

def pose_from_amber_params(crd_path, top_path, index=0):
    '''Takes an AMBER coordinates and topology file and creates a pose from them
    with `TrajToPoses`_. If the coordinates contain multiple frames, it takes
    either the first one, or the one specified by the argument ``index``.

    Parameters
    ----------
    crd_path : str
        Path to the input .rst7 or .nc file.
    top_path : str
        Path to the input .parm7 file.
    index : int, optional
        The index of the frame to use from the trajectory. 0 by default.

    Returns
    -------
    rosetta.core.Pose
        The pose corresponding to the given input files.
    '''

    return TrajToPoses(pt.iterload(crd_path, top_path))[index]
