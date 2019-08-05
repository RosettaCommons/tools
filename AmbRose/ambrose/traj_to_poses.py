'''The `TrajToPoses`_ class and its supporting functions.'''

# pragma pylint: disable=bad-whitespace

import collections, enum
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

class _ResidueClass(enum.Flag):
    '''Describes attributes of an AMBER residue, used to control state
    transition for the _TopologyParser finite state machine. The possible flags
    are:

    TER
        Special flag reserved for the end of a topology. The only flag that
        doesn't represent a residue.
    NAA
        Amino acid at the N-terminal of a protein.
    AA
        Amino acid at neither terminal of a protein.
    CAA
        Amino acid at the C-terminal of a protein.
    RNA5
        RNA base at the 5' end of an RNA strand.
    RNA
        RNA base at neither end of an RNA strand.
    RNA3
        RNA base at the 3' end of an RNA strand.
    DNA5
        DNA base at the 5' end of an RNA strand.
    DNA
        DNA base at neither end of an RNA strand.
    DNA3
        DNA base at the 3' end of an RNA strand.
    ACE
        An N-terminal acetyl on a protein.
    NME
        A C-terminal N-methylamine on a protein.

    Additionally, the following composite flags (flags that are combinations of
    single flags) are defined:

    NONE
        Something that doesn't belong to any class.
    NAA_OR_CAA
        Amino acid at either terminal of a protein.
    ANY_AA
        Amino acid.
    ANY_RNA
        RNA base.
    ANY_DNA
        DNA base.
    NA5
        RNA or DNA base at the 5' end of a strand.
    NA
        RNA or DNA base at neither end of a strand.
    NA3
        RNA or DNA base at the 5' end of a strand.
    NA5_OR_NA3
        RNA or DNA base at either end of a strand.
    ANY_NA
        Nucleobase.
    ROSETTA_EQUIVALENT
        Any residue that has a full Rosetta residue named after it; currently
        this is amino acids and nucleobases.
    '''
    NONE       = 0
    TER        = enum.auto()
    NAA        = enum.auto()
    AA         = enum.auto()
    CAA        = enum.auto()
    NAA_OR_CAA = NAA | CAA
    ANY_AA     = NAA | AA | CAA
    RNA5       = enum.auto()
    RNA        = enum.auto()
    RNA3       = enum.auto()
    ANY_RNA    = RNA5 | RNA | RNA3
    DNA5       = enum.auto()
    DNA        = enum.auto()
    DNA3       = enum.auto()
    ANY_DNA    = DNA5 | DNA | DNA3
    NA5        = RNA5 | DNA5
    NA         = RNA | DNA
    NA3        = RNA3 | DNA3
    NA5_OR_NA3 = NA5 | NA3
    ANY_NA     = NA5 | NA | NA3
    ACE        = enum.auto()
    NME        = enum.auto()
    ROSETTA_EQUIVALENT = ANY_AA | ANY_NA

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
    name : str
        The 4-character name of the atom in its Rosetta ``ResidueType``.
    res_name : str
        The 3-4-character name of the residue the atom belongs to.
    chain_id : str
        The 1-character ID of the chain the atom belongs to.
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

    def __init__(self, name, res_name, chain_id, element, basep=True):
        self.name     = name
        self.res_name = res_name
        self.chain_id = chain_id
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

        ai = pr.rosetta.core.io.AtomInformation()
        ai.resName   = self.res_name
        ai.name      = self.name
        ai.chainID   = self.chain_id
        ai.occupancy = 1.0
        ai.segmentID = '    '
        ai.element   = self.element
        return ai

class _AMBERResidue:
    '''Represents an AMBER residue, with lots of derived information about it
    that's relevant to mapping it to a Rosetta residue.

    This class emphatically does *not* store the original ``Residue`` object,
    and its four listed attributes are all mutable. This is to provide
    flexibility in the case that weird hacks have to be done to make edge cases
    work later. (I have only a small but measurable amount of shame about this.)

    Parameters
    ----------
    top : pytraj.Topology
        The topology to which the residue belongs.
    residue_i : int
        The index of the Pytraj residue object that we're casting within its
        topology.

    Attributes
    ----------
    i : int
        The index of the residue.
    name : str
        The full name of the residue.
    regular_name : str
        The full name of the residue, minus the N or C at the front for
        N/C-terminal amino acids, and minus the 5 or 3 at the end for 5'/3'
        nucleic acids. Similar to the representative type name for Rosetta
        residues.
    first : int
        The index of the first atom in the residue.
    last : int
        The index of the last atom in the residue.
    atoms : tuple[str]
        The names of the residue's atoms.
    res_class : `_ResidueClass`_
        The residue class of the residue.
    '''

    # constant sets of residue names for the classifier
    RES_CLASS_MAP = {
        frozenset({'ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'CYM', 'CYS', 'CYX',
                   'GLN', 'GLU', 'GLH', 'GLY', 'HID', 'HIE', 'HIP', 'HYP',
                   'ILE', 'LEU', 'LYN', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                   'THR', 'TRP', 'TYR', 'VAL'}):
            _ResidueClass.AA,
        frozenset({'A', 'C', 'G', 'U'}): _ResidueClass.RNA,
        frozenset({'DA', 'DC', 'DG', 'DT'}): _ResidueClass.DNA,
    }
    # Maps of standard AMBER residue names to Rosetta residue names, with the
    # caveat that the names with N or C prepended for N-terminal or C-terminal
    # residues are not included. Only entries where the name is actually
    # different are included.
    RES_NAME_MAP = {
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
        'C'  :'RCY',
        'G'  :'RGU',
        'U'  :'URA',
        'DA' :'ADE',
        'DC' :'CYT',
        'DG' :'GUA',
        'DT' :'THY',
    }

    def __init__(self, top, res_i):
        self.i = res_i
        residue = top.residue(res_i)
        self.name = residue.name
        self.first = residue.first
        self.last = residue.last
        self.atoms = tuple(a.name for a in top[self.first:self.last].atoms)
        self.res_class = None
        self.regular_name = self.name
        self.patches = []
        for name_set, res_class in _AMBERResidue.RES_CLASS_MAP.items():
            if self.name in name_set:
                self.res_class = res_class
                break
            if res_class & _ResidueClass.AA \
               and len(self.name) == 4 \
               and self.name[-3:] in name_set \
               and self.name[0] in 'NC':
                self.regular_name = self.name[-3:]
                self.res_class = res_class
                if self.name[0] == 'N':
                    self.is_naa = True
                elif self.name[0] == 'C':
                    self.is_caa = True
                break
            if res_class & _ResidueClass.NA \
               and len(self.name) in (2, 3) \
               and self.name[:-1] in name_set \
               and self.name[-1] in '53':
                self.regular_name = self.name[:-1]
                self.res_class = res_class
                if self.name[-1] == '5':
                    self.is_na5 = True
                elif self.name[-1] == '3':
                    self.is_na3 = True
                break
        if self.is_aa:
            if 'OXT' in self.atoms:
                self.is_aa = False
                self.is_caa = True
            else:
                nitro_index = self.first + self.atoms.index('N')
                nitro_bonded_names = \
                    tuple(top.atom(
                              tuple(b.indices)[
                                  not tuple(b.indices).index(nitro_index)
                               ]).name \
                          for b in top.bonds if nitro_index in b.indices)
                if 'C' not in nitro_bonded_names:
                    self.is_aa = False
                    self.is_naa = True
        self.rosetta_name = _AMBERResidue.RES_NAME_MAP.get(
            self.regular_name,
            self.regular_name) if self.has_rosetta_equivalent else None
    def __repr__(self):
        return '{} {} {}'.format(self.name, self.i, self.res_class)
    @property
    def n_atoms(self):
        '''int: Number of atoms in residue.'''
        return self.last - self.first + 1
    def has_class(self, res_class):
        return self.res_class & res_class
    def add_class_and_patch(self, res_class=_ResidueClass.NONE, patch=None):
        self.res_class |= res_class
        if patch is not None:
            self.patches.append(patch)
    def remove_class_and_patch(self, res_class=_ResidueClass.NONE, patch=None):
        self.res_class &= ~res_class
        if patch is not None:
            try:
                self.patches.remove(patch)
            except ValueError:
                pass
    @property
    def is_aa(self):
        return self.has_class(_ResidueClass.AA)
    @is_aa.setter
    def is_aa(self, value):
        CLASS = _ResidueClass.AA
        if value:
            self.add_class_and_patch(CLASS)
        else:
            self.remove_class_and_patch(CLASS)
    @property
    def is_naa(self):
        return self.has_class(_ResidueClass.NAA)
    @is_naa.setter
    def is_naa(self, value):
        CLASS = _ResidueClass.NAA
        PATCH = 'NtermProteinFull'
        if value:
            self.add_class_and_patch(CLASS, PATCH)
        else:
            self.remove_class_and_patch(CLASS, PATCH)
    @property
    def is_caa(self):
        return self.has_class(_ResidueClass.CAA)
    @is_caa.setter
    def is_caa(self, value):
        CLASS = _ResidueClass.CAA
        PATCH = 'CtermProteinFull'
        if value:
            self.add_class_and_patch(CLASS, PATCH)
        else:
            self.remove_class_and_patch(CLASS, PATCH)
    @property
    def is_naa_or_caa(self):
        return self.has_class(  _ResidueClass.NAA \
                              | _ResidueClass.CAA)
    @property
    def is_na(self):
        return self.has_class(_ResidueClass.NA)
    @is_na.setter
    def is_na(self, value):
        CLASS = _ResidueClass.NA
        if value:
            self.add_class_and_patch(CLASS)
        else:
            self.remove_class_and_patch(CLASS)
    @property
    def is_na5(self):
        return self.has_class(_ResidueClass.NA5)
    @is_na5.setter
    def is_na5(self, value):
        CLASS = _ResidueClass.NA5
        PATCH = 'LowerRNA'
        if value:
            self.add_class_and_patch(CLASS, PATCH)
        else:
            self.remove_class_and_patch(CLASS, PATCH)
    @property
    def is_na3(self):
        return self.has_class(_ResidueClass.NA3)
    @is_na3.setter
    def is_na3(self, value):
        CLASS = _ResidueClass.NA3
        PATCH = 'UpperRNA'
        if value:
            self.add_class_and_patch(CLASS, PATCH)
        else:
            self.remove_class_and_patch(CLASS, PATCH)
    @property
    def is_na5_or_na3(self):
        return self.has_class(  _ResidueClass.NA5 \
                              | _ResidueClass.NA3)
    @property
    def has_rosetta_equivalent(self):
        return self.has_class(_ResidueClass.ROSETTA_EQUIVALENT)

class _RosettaResidue:
    '''Describes the key properties of a ResidueType used to construct
    :obj:`_AtomRecord`s, and a mapping of its atoms onto absolute indices in a
    pytraj ``Topology``. This object can be indexed and iterated over to
    retreive its :obj:`_AtomRecord`s, but not sliced.

    Parameters
    ----------
    amber_residue : `_AMBERResidue`_
        The AMBER residue it corresponds to. While some residues (such as
        acetylated N-terminal residues) may correspond to multiple residues in
        AMBER, only one of them will have the property ``.is_rosetta_res``. This
        is generally the biggest, most "canonical" residue.
    chain_id : str, optional
        A single printable ASCII character representing the chain to which this
        residue belongs. '@' by default (0x40).

    Attributes
    ----------
    amber_residue : `_AMBERResidue`_
        The AMBER residue representing this residue.
    chain_id : str
        A single printable ASCII character representing the chain to which this
        residue belongs.
    atoms : list(:obj:`_AtomRecord`)
        A list of the :obj:`_AtomRecord`s of this residue.
    __names_to_indices : dict
        A dict that maps the 4-character names of the atoms in the corresponding
        Rosetta ``ResidueType`` to their indices in `atoms`.
    '''

    def __init__(self, amber_residue, chain_id='@'):
        name = amber_residue.rosetta_name
        patches = amber_residue.patches
        res_type = _get_pr_res(':'.join((name, *patches)))

        ## create _AtomRecords for self.atoms
        self.atoms = []
        atom_names_to_indices = {}
        all_base_atoms =    amber_residue.is_naa_or_caa \
                         or amber_residue.is_na5_or_na3
        for i in range(1, res_type.natoms()):
            atom_name = res_type.atom_name(i)
            atom_names_to_indices[atom_name] = i
            self.atoms.append(
                _AtomRecord(atom_name,
                            res_type.name3(),
                            chain_id,
                            str(res_type.atom(i).element())[9:],
                            basep=(all_base_atoms or \
                                   res_type.has(atom_name))))
        self.atoms = tuple(self.atoms)

        def set_indices_by_name_relative(indices_dict):
            '''Provided a dict that maps Rosetta atom names to an offset of its
            pt_i from the index of the first atom in the AMBER residue, set the
            pt_is of each named atom to the appropriate absolute index.

            Parameters
            ----------
            indices_dict : dict
                A dict mapping the Rosetta names of atoms in this residue to the
                indices of their corresponding atoms in a pytraj trajectory,
                relative to the beginning of the base AMBER residue.
            '''

            for atom_name, i in indices_dict.items():
                self[atom_names_to_indices[atom_name]].pt_i = \
                    amber_residue.first + i

        ## assign indices for N-terminal acetyl atoms
        if 'N_acetylated' in patches:
            set_indices_by_name_relative(
                {' CP ': -2,
                 ' CQ ': -5,
                 ' OCP': -1})

        ## assign indices for hydrogens of base residue
        # This sets the pt_i of the base atoms, based on a map of AMBER atom
        # names to absolute pt_is, via a specific algorithm that works for
        # canonical amino acids and canonical nucleic acids (and
        # hydroxyproline, and not much else).
        pt_names_to_indices = {atom_name: amber_residue.first + i \
                               for i, atom_name \
                               in enumerate(amber_residue.atoms)}
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
                if name in ('PRO', 'HPR') and atom_record.name == ' NV ':
                    atom_record.h_number = None
                    atom_record.coda = 'N'

                ## Actually retrieve index:
                atom_record.pt_i = pt_names_to_indices.get(atom_record.coda)
                # we actually have no plan if the .get returns None
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

        ## assign indices for C-terminal N-methylamine atoms
        if 'C_methylamidated' in patches:
            n = amber_residue.n_atoms()
            set_indices_by_name_relative(
                {' NR ': n,
                 ' CS ': n + 2,
                 ' HR ': n + 1})

    def __len__(self):
        return len(self.atoms)
    def __iter__(self):
        for i in range(len(self)):
            yield self[i]
    def __getitem__(self, key):
        return self.atoms[key]

class _TopologyParser:
    '''A finite state machine that parses residues from a pytraj topology and
    builds a list of `_RosettaResidue`_ objects.'''

    def __init__(self):
        ## Build transitions:
        # Each transition emits a tuple of "action tokens" that cause the parser
        # to perform some action. Note that a "Rosetta equivalent queue" exists,
        # exists, to which every AMBER residue in the topology that has a
        # Rosetta equivalent immediately gets pushed, before any actions are
        # performed. The tokens are:
        #
        # SKIP
        #     If the current residue is a Rosetta equivalent residue, remove it
        #     from the queue. Otherwise do nothing.
        # CHBR
        #     Increment the chain counter. Whenever a residue is EMITted, it
        #     receives whatever chain ID the current chain counter reads.
        # EMIT
        #     Package the first residue in the queue into a `_RosettaResidue`_
        #     object and add it to the list of those objects. This removes the
        #     residue from the queue, and no further actions may be performed
        #     on it. (Unless in the future an action token is invented that
        #     operates on the Rosetta residue list, but let's not.)
        # PATCH:patch_name
        #     Adds the patch *patch_name* to the end of the list of patches to
        #     be applied to the residue at the front of the queue. *patch_name*
        #     may be any string, but in practice should be the name of a valid
        #     Rosetta patch.
        self.state = 'BEGIN'
        Class = _ResidueClass # brevity
        self.TRANSITIONS = \
            {'BEGIN':  {Class.NAA: ('PRE_AA', ('CHBR', 'EMIT',)),
                        Class.ACE: ('ACE_AA', ('CHBR',)),
                        Class.NA5: ('NA',     ('CHBR', 'EMIT',))},
             'PRE_AA': {Class.AA:  ('AA',     ()),
                        Class.CAA: ('BEGIN',  ('EMIT',))},
             'ACE_AA': {Class.AA:  ('AA',     ('PATCH:N_acetylated',))},
             'AA':     {Class.AA:  ('AA',     ('EMIT',)),
                        Class.CAA: ('BEGIN',  ('EMIT', 'EMIT',)),
                        Class.NME: ('BEGIN',  ('PATCH:C_methylamidated', 'EMIT'))},
             'NA':     {Class.NA:  ('NA',     ('EMIT',)),
                        Class.NA3: ('BEGIN',  ('EMIT',))}}
        self.DEFAULT_TRANSITIONS = \
            {'BEGIN':  ('BEGIN', ('SKIP',)),
             'ACE_AA': ('AA',    ('SKIP',)),
             'AA':     ('AA',    ('SKIP',)),
             'NA':     ('NA',    ('SKIP',))}
    def transit(self, received_class):
        '''Perform a state transition, with received_id as the key that maps to
        the transition.'''

        for accepted_class, possible_transition \
            in self.TRANSITIONS[self.state].items():
            if received_class & accepted_class:
                next_state, emitted = possible_transition
                break
        else: # if we went through the whole loop without breaking...
            next_state, emitted = self.DEFAULT_TRANSITIONS[self.state]
        self.state = next_state
        return emitted
    def parse_topology(self, top):
        '''Returns a list of `_ResidueMap`_ objects, which contain the
        information needed to build ``core::chemical::Residue`` objects from the
        list of atoms in the given topology, namely ResidueType names, any
        Patches, and a map of each ResidueType's atoms to an atom index in the
        given topology.'''

        amber_res = None
        output_residues = []          # collects output _RosettaResidues
        rosetta_equivalent_queue = [] # queues Rosetta-equivalent _AMBERResidues
        chain_id = '@'

        def process_tokens(tokens):
            nonlocal amber_res
            nonlocal output_residues
            nonlocal rosetta_equivalent_queue
            nonlocal chain_id
            for token in tokens:
                utils.debug_print('token:', token)
                if token == 'SKIP':
                    if amber_res is not None \
                       and amber_res.has_rosetta_equivalent:
                        del rosetta_equivalent_queue[-1]
                    continue
                if token == 'CHBR':
                    if ord(chain_id) < 126:
                        chain_id = chr(ord(chain_id) + 1)
                    continue
                if token.startswith('PATCH:'):
                    try:
                        rosetta_equivalent_queue[0].add_class_and_patch(
                            patch=token[6:])
                    except IndexError:
                        raise errors.TopologyParsingError(
                            'Error applying patch {} to topology '
                            'at residue: {}'.format(token, amber_res))
                    continue
                if token == 'EMIT':
                    try:
                        output_residues.append(
                            _RosettaResidue(
                                rosetta_equivalent_queue.pop(0),
                                chain_id))
                    except IndexError:
                        raise errors.TopologyParsingError(
                            'Error applying adding residue to topology '
                            'at residue: {}'.format(amber_res))
                    continue

        for res_i in range(top.n_residues):
            amber_res = _AMBERResidue(top, res_i)
            utils.debug_print('residue:', amber_res)
            if amber_res.has_rosetta_equivalent:
                rosetta_equivalent_queue.append(amber_res)
            process_tokens(self.transit(amber_res.res_class))

        amber_res = None
        process_tokens(self.transit(_ResidueClass.TER))
        return tuple(output_residues)

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
    __residues : tuple(_ResidueMap), optional
        A tuple of `_RosettaResidue`_s to use as its internal map of Rosetta
        residues to pytraj residues. Used whenever slices of this object are
        made.

    Attributes
    ----------
    traj : pytraj.Trajectory
        The trajectory the Poses are derived from.
    options : rosetta.core.import_pose.ImportPoseOptions
        The options used to derive the Poses.
    __residues : tuple(_ResidueMap)
        A tuple of `_RosettaResidue`_s to use as its internal map of Rosetta
        residues to pytraj residues.
    '''

    def __init__(self, traj, options=None, *, __residues=None):
        self.traj = traj
        self.options = options or pr.rosetta.core.import_pose.ImportPoseOptions
        if __residues is not None:
            self.__residues = __residues
        else:
            self.__residues = _TopologyParser().parse_topology(traj.top)
    def __len__(self):
        return len(self.traj)
    def __iter__(self):
        for i in range(len(self)):
            yield self[i]
    def __getitem__(self, key):
        if isinstance(key, slice):
            return TrajToPoses(self.traj[key], __residues=self.__residues)

        frame = self.traj[key]
        sfr  = pr.rosetta.core.io.StructFileRep()
        ChainAtoms      = pr.rosetta.utility.vector0_core_io_AtomInformation
        chain_atoms_map = collections.defaultdict(ChainAtoms)
        atom_i = itertools.count(1)
        for res_i, res in enumerate(self.__residues, 1):
            for atom_record in res:
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
