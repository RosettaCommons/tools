import sys
import PSSM
from Bio.PDB import *
from Bio import pairwise2
from rosettautil.util import fileutil
import math
import warnings
from Bio import SVDSuperimposer
from itertools import izip
import amino_acids

def sequence_recovery(native_struct, designed_struct):
    """calculate percent sequence recovery between a native and designed struct"""
    native_residues = native_struct.get_residues()
    designed_residues = designed_struct.get_residues()
    total = 0.0
    recovered = 0.0
    for native, designed in zip(native_residues, designed_residues):
        if native.get_resname() == designed.get_resname():
            recovered += 1
        total += 1
    return recovered / total


def sequence_recovery_allowed_list(native_struct, designed_struct, allowed_residues):
    """calculate percent sequence recovery between a native and designed struct"""
    native_residues = native_struct.get_residues()
    designed_residues = designed_struct.get_residues()
    total = 0.0
    recovered = 0.0
    for native, designed in zip(native_residues, designed_residues):
        full_id = native.get_full_id()
        residue_data = (full_id[3][1], full_id[2])
        if residue_data not in allowed_residues:
            continue
        if native.get_resname() == designed.get_resname():
            recovered += 1
        total += 1
    return recovered / total


def sequence_recovery_range(native_struct, designed_struct, min, max):
    """calculate percent sequence recovery between a native and designed
    struct for all residues between a minimum and maximum b factor """
    native_residues = native_struct.get_residues()
    designed_residues = designed_struct.get_residues()
    total = 0.0
    recovered = 0.0
    for native, designed in zip(native_residues, designed_residues):
        score = native.get_list()[1].get_bfactor()
        if score >= min and score <= max:
            if native.get_resname() == designed.get_resname():
                recovered += 1
            total += 1
    return recovered / total


def sequence_recovery_group(native_struct, designed_struct, min, max):
    """calculate sequence recovery by group (non-polar, polar,
    aromatic, and charged) for all residues between a minimum and
    maximum b factor"""
    non_polar = ["GLY", "ALA", "VAL", "LEU", "MET", "ILE"]
    polar = ["SER", "THR", "CYS", "PRO", "ASN", "GLN"]
    aromatic = ["PHE", "TYR", "TRP"]
    charged = ["LYS", "ARG", "HIS", "ASP", "GLU"]

    groups = {"non_polar": non_polar, "polar": polar, "aromatic": aromatic, "charged": charged}

    native_residues = native_struct.get_residues()
    designed_residues = designed_struct.get_residues()
    total = {"non_polar": 0.0, "polar": 0.0, "aromatic": 0.0, "charged": 0.0}
    recovered = {"non_polar": 0.0, "polar": 0.0, "aromatic": 0.0, "charged": 0.0}

    for native, designed in zip(native_residues, designed_residues):
        score = native.get_list()[1].get_bfactor()

        if score >= min and score <= max:
            for group in groups:
                if native.get_resname() in groups[group] and designed.get_resname() in groups[group]:
                    recovered[group] += 1
                if designed.get_resname() in groups[group]:
                    total[group] += 1
                    break
    for group in groups:
        print group, recovered[group], total[group]
        # recovered[group]=recovered[group]/total[group]

    return recovered


def sequence_composition_allowed_list(struct, allowed_residues):
    """calculate sequence composition by residue"""
    struct_residues = struct.get_residues()
    composition = {}
    for residue in struct_residues:
        full_id = residue.get_full_id()
        residue_data = (full_id[3][1], full_id[2])
        if residue_data not in allowed_residues:
            continue
        residue_name = residue.get_resname()
        try:
            composition[residue_name] += 1
        except KeyError:
            composition[residue_name] = 1
    return composition


def sequence_composition(struct):
    """calculate sequence composition by residue"""
    struct_residues = struct.get_residues()
    composition = {}
    for residue in struct_residues:
        residue_name = residue.get_resname()
        try:
            composition[residue_name] += 1
        except KeyError:
            composition[residue_name] = 1
    return composition


def sequence_composition_range(struct, min, max):
    """calculate sequence composition for all residues within a
    minimum and maximum b factor"""
    struct_residues = struct.get_residues()
    composition = {}
    for residue in struct_residues:
        score = residue.get_list()[1].get_bfactor()
        if score >= min and score <= max:
            residue_name = residue.get_resname()
            try:
                composition[residue_name] += 1
            except KeyError:
                composition[residue_name] = 1
    return composition


def pssm_recovery_map_allowed_list(native_struct, designed_struct, pssm_map, allowed_residues):
    """calculate the pssm recovery given a structure and a pssm map"""
    native_residues = native_struct.get_residues()
    designed_residues = designed_struct.get_residues()
    # pssm_recovery = 0.0;
    # struct_size = 0.0;
    recovery_map = {}
    for native, designed in zip(native_residues, designed_residues):

        full_id = native.get_full_id()
        residue_data = (full_id[3][1], full_id[2])
        if residue_data not in allowed_residues:
            continue
        designed_name = designed.get_resname()
        native_name = native.get_resname()
        designed_num = designed.get_id()[1]
        try:
            status = pssm_map.conserved(designed_num, designed_name)
        except KeyError:
            warnings.warn("ignoring noncanonical amino acid " + residue_name + " in pssm calculation")
            continue
        if status:
            try:
                recovery_map[native_name] += 1
            except KeyError:
                recovery_map[native_name] = 1
    return recovery_map


def pssm_recovery_map(struct, pssm_map):
    """calculate the pssm recovery given a structure and a pssm map"""
    struct_residues = struct.get_residues()
    # pssm_recovery = 0.0;
    # struct_size = 0.0;
    recovery_map = {}
    for residue in struct_residues:
        residue_name = residue.get_resname()
        residue_num = residue.get_id()[1]
        try:
            status = pssm_map.conserved(residue_num, residue_name)
        except KeyError:
            warnings.warn("ignoring noncanonical amino acid " + residue_name + " in pssm calculation")
            continue
        if status:
            try:
                recovery_map[residue_name] += 1
            except KeyError:
                recovery_map[residue_name] = 1
    return recovery_map


def pssm_recovery_map_range(struct, pssm_map, min, max):
    """calculate the pssm recovery within a range of b factors given a
    structure and a pssm map"""
    struct_residues = struct.get_residues()
    recovery_map = {}
    for residue in struct_residues:
        score = residue.get_list()[1].get_bfactor()
        if score >= min and score <= max:
            residue_name = residue.get_resname()
            residue_num = residue.get_id()[1]
            try:
                status = pssm_map.conserved(residue_num, residue_name)
            except KeyError:
                warnings.warn("ignoring noncanonical amino acid " + residue_name + " in pssm calculation")
                continue
            if status:
                try:
                    recovery_map[residue_name] += 1
                except KeyError:
                    recovery_map[residue_name] = 1
    return recovery_map


def pssm_recovery_allowed_list(designed_struct, allowed_residues, pssm_map, manual_size=0):
    """return percent pssm recovery"""
    designed_residues = designed_struct.get_residues()
    pssm_recovery = 0.0
    struct_size = 0.0
    for residue in designed_residues:

        full_id = residue.get_full_id()
        residue_data = (full_id[3][1], full_id[2])
        if residue_data not in allowed_residues:
            continue
        residue_name = residue.get_resname()
        residue_num = residue.get_id()[1]
        try:
            status = pssm_map.conserved(residue_num, residue_name)
        except KeyError:
            warnings.warn("ignoring noncanonical amino acid " + residue_name + " in pssm calculation")
            continue
        if status:
            pssm_recovery += 1.0
        struct_size += 1.0
    if manual_size == 0:
        return pssm_recovery / struct_size
    else:
        return pssm_recovery / float(manual_size)


def pssm_recovery(struct, pssm_map, manual_size=0):
    """return percent pssm recovery"""
    struct_residues = struct.get_residues()
    pssm_recovery = 0.0
    struct_size = 0.0
    for residue in struct_residues:
        residue_name = residue.get_resname()
        residue_num = residue.get_id()[1]
        try:
            status = pssm_map.conserved(residue_num, residue_name)
        except KeyError:
            warnings.warn("ignoring noncanonical amino acid " + residue_name + " in pssm calculation")
            continue
        if status:
            pssm_recovery += 1.0
        struct_size += 1.0
    if manual_size == 0:
        return pssm_recovery / struct_size
    else:
        return pssm_recovery / float(manual_size)


def pssm_recovery_range(struct, pssm_map, min, max):
    """return percent pssm recovery fro residues within a range of b factors"""
    pssm_recovery = 0.0
    struct_size = 0.0
    for residue in struct.get_residues():
        score = residue.get_list()[1].get_bfactor()
        # print score
        if score >= min and score <= max:
            residue_name = residue.get_resname()
            residue_num = residue.get_id()[1]
            try:
                status = pssm_map.conserved(residue_num, residue_name)
            except KeyError:
                warnings.warn("ignoring noncanonical amino acid " + residue_name + " in pssm calculation")
                continue
            if status:
                pssm_recovery += 1.0
            struct_size += 1.0

    return pssm_recovery / struct_size


def pssm_scores(struct, pssm_map):
    """ return raw pssm total for each residue"""
    struct_residues = struct.get_residues()
    pssm_scores = {}
    size = 0
    for residue in struct_residues:
        size += 1
        residue_name = residue.get_resname()
        residue_num = residue.get_id()[1]
        try:
            pssm_scores[residue_name] += pssm_map.get_score(residue_num, residue_name)
        except KeyError:
            pssm_scores[residue_name] = pssm_map.get_score(residue_num, residue_name)
    for key in pssm_scores:
        pssm_scores[key] = pssm_scores[key] / size
    return pssm_scores


def pssm_scores_range(struct, pssm_map, min, max):
    """ return raw pssm total for each residue within a range of pssm scores"""
    struct_residues = struct.get_residues()
    pssm_scores = {}
    size = 0
    for residue in struct_residues:
        score = residue.get_list()[1].get_bfactor()
        # print score
        if score >= min and score <= max:
            size += 1
            residue_name = residue.get_resname()
            residue_num = residue.get_id()[1]
            try:
                pssm_scores[residue_name] += pssm_map.get_score(residue_num, residue_name)
            except KeyError:
                pssm_scores[residue_name] = pssm_map.get_score(residue_num, residue_name)
    for key in pssm_scores:
        pssm_scores[key] = pssm_scores[key] / size
    return pssm_scores


def ca_rms_only(struct_a, struct_b, residue_list):
    """calculate the CA rmsd of two structures using the residues in residue_list"""
    residues_a = struct_a.get_residues()
    residues_b = struct_b.get_residues()

    d_2_sum = 0.0
    resn = 0
    for (res_a, res_b) in zip(residues_a, residues_b):
        if res_a.get_id()[1] not in residue_list:
            continue
        CA_a = res_a['CA']
        CA_b = res_b['CA']

        distance_2 = (CA_a - CA_b) ** 2
        d_2_sum += distance_2
        resn += 1

    rmsd = math.sqrt(d_2_sum / resn)
    return rmsd


def atom_rms(atoms_a, atoms_b, debug=False, force=False):
    """calculate the all atom rmsd given two lists of atoms and a list of residues"""
    garbage_counter = 0
    included_counter = 0
    d_2_sum = 0.0
    resn = 0
    if debug:
        debug_set_a = []
        debug_set_b = []
    for(atom_a, atom_b) in zip(atoms_a, atoms_b):
        # print "calculating for", parent_a.get_id()
        distance_2 = (atom_a - atom_b) ** 2
        if atom_a.get_id() != atom_b.get_id():
            if force is True:
                print str(atom_a.id) + "_" + str(atom_a.get_parent().id[1]) + "_" + str(atom_a.get_parent().get_parent().id), \
                    "and", str(atom_b.id) + "_" + str(atom_b.get_parent().id[1]) + "_" + str(atom_b.get_parent().get_parent(
                    ).id), "don't match up, tossing them out of the rmsd calculation"
                garbage_counter += 1
                continue
            else:
                raise Exception("Atom ordering is not the same between the two pdbs, atom ids don't match")
        else:
            included_counter += 1
        if debug:
            debug_set_a.append(atom_a)
            debug_set_b.append(atom_b)
        d_2_sum += distance_2
        resn += 1
    if debug:
        print "removed %i atom sets from RMSD calculation" % (garbage_counter)
        print "calculating rmsd on %i sets of atoms" % (included_counter)
    rmsd = math.sqrt(d_2_sum / resn)
    if debug:
        with open("debug_atoms.txt", 'w') as f:
            for i, j in zip(debug_set_a, debug_set_b):
                atom_i, atom_j = i.get_id(), j.get_id()
                res_i, res_j = i.get_parent().get_id()[1], j.get_parent().get_id()[1]
                chain_i, chain_j = i.get_parent().get_parent().get_id(), j.get_parent().get_parent().get_id()
                f.write("chain {0},residue {1},atom {2} v chain {3}, residue {4}, atom {5}\n".format(
                    chain_i, res_i, atom_i, chain_j, res_j, atom_j))
    return rmsd


def copy_b_factor(native_pdb, designed_pdb):
    """ copy the b factors from one structure to another"""
    native_atoms = native_pdb.get_atoms()
    designed_atoms = designed_pdb.get_atoms()
    for native_atom, designed_atom in zip(native_atoms, designed_atoms):
        designed_atom.set_bfactor(native_atom.get_bfactor)
    return designed_pdb


def calculate_rms(native, decoy, ca_mode, residues, rms_residues, chain, debug, force):
    native_atoms = []
    decoy_atoms = []
    residue_set = set()
    rms_residue_set = set()
    # check to see if we are calculating rmsd on just residues
    if residues:
        residue_file = open(residues, 'r')
        for residue in residue_file:
            residue = residue.strip()
            if residue is not '':
                residue_set.add(int(residue))
        residue_file.close()

    # this seems redundant
    if rms_residues:
        rms_residue_file = open(rms_residues, 'r')
        for rms_residue in rms_residue_file:
            rms_residue = rms_residue.strip()
            if rms_residue is not '':
                rms_residue_set.add(int(rms_residue))
            rms_residue_file.close()

    # only get atoms of a certain chain
    for(native_chain, decoy_chain) in zip(native[0].get_list(), decoy[0].get_list()):
        if chain:
            if native_chain.get_id() not in chain:
                print "ignoring chain " + native_chain.get_id()
                continue
    # iterate through residues
        for(native_residue, decoy_residue) in zip(native_chain.get_list(), decoy_chain.get_list()):
                # check to see if it inclues our residues
            if len(residue_set) > 0 and native_residue.id[1] in residue_set:
                if ca_mode == "ca":
                    try:
                        native_atoms.append(native_residue['CA'])
                        decoy_atoms.append(decoy_residue['CA'])
                    except KeyError:
                        print "WARNING: residue", str(native_residue.get_id()[1]), "has no CA atom in either the native or the decoy structure.  Either this residue is a hetatm or part of your backbone is missing. The residue will not be included in the RMSD calculations."
                elif ca_mode == "bb":
                    try:
                        native_atoms.append(native_residue['CA'])
                        native_atoms.append(native_residue['C'])
                        native_atoms.append(native_residue['O'])
                        native_atoms.append(native_residue['N'])
                        decoy_atoms.append(decoy_residue['CA'])
                        decoy_atoms.append(decoy_residue['C'])
                        decoy_atoms.append(decoy_residue['O'])
                        decoy_atoms.append(decoy_residue['N'])
                    except KeyError:
                        print "WARNING: residue", str(native_residue.get_id()[1]), "has no backbone atoms"
                elif ca_mode == "all":
                    for (native_atom, decoy_atom) in zip(native_residue.get_list(), decoy_residue.get_list()):
                        native_atoms.append(native_atom)
                        decoy_atoms.append(decoy_atom)
                else:
                    print "you must specify a ca_mode, all,bb, or ca"
                    sys.exit()

            # no residue set
            elif not residue_set:
                if ca_mode == 'ca':
                    try:
                        native_atoms.append(native_residue['CA'])
                        decoy_atoms.append(decoy_residue['CA'])
                    except KeyError:
                        print "WARNING: residue", str(native_residue.get_id()[1]), "has no CA atom in either the native or the decoy structure.  Either this residue is a hetatm or part of your backbone is missing. The residue will not be included in the RMSD calculations."
                elif ca_mode == 'bb':
                    try:
                        native_atoms.append(native_residue['CA'])
                        native_atoms.append(native_residue['C'])
                        native_atoms.append(native_residue['N'])
                        native_atoms.append(native_residue['O'])
                        decoy_atoms.append(decoy_residue['CA'])
                        decoy_atoms.append(decoy_residue['C'])
                        decoy_atoms.append(decoy_residue['N'])
                        decoy_atoms.append(decoy_residue['O'])
                    except KeyError:
                        print "WARNING: residue", str(native_residue.get_id()[1]), "has no backbone atoms"
                elif ca_mode == 'all':
                    for (native_atom, decoy_atom) in zip(native_residue.get_list(), decoy_residue.get_list()):
                        native_atoms.append(native_atom)
                        decoy_atoms.append(decoy_atom)
                else:
                    print "you must specify a ca_mode, all,bb, or ca"
                    sys.exit()
    superpose = Superimposer()
    superpose.set_atoms(native_atoms, decoy_atoms)
    superpose.apply(decoy.get_atoms())
    return atom_rms(native_atoms, decoy_atoms, debug, force)


def calculate_all_superpositions(native, decoy, residues, debug,pairwise):
    """calculates c alpha , backbone and all atom rmsds. you can pass a list of residues if you want. pass the resfile
    through the resfiles class first which returns a resfile list:\n
    1 g allaa\n
    3 h pikaa g\n
    would return a list of tuples like :\n
    [(1,g),(3,h)]"""
    native_atoms_ca = []
    decoy_atoms_ca = []
    native_atoms_bb = []
    decoy_atoms_bb = []
    native_atoms_all = []
    decoy_atoms_all = []
    # check to see if we are calculating rmsd on just residues
    counter = 1
    what_to_skip_native = []
    what_to_skip_decoy = []
    native_residues = [i for i in native.get_residues()]
    decoy_residues = [i for i in decoy.get_residues()]

    #first figure out what to skip
    for fasta1,fasta2 in zip(pairwise['native'],pairwise['target']):
            if fasta1 == "-":
                    what_to_skip_native.append(counter)
            elif fasta2 == "-":
                    what_to_skip_decoy.append(counter)
            counter += 1

    for counter,residue in enumerate(native.get_residues(),start=1):
        if counter in what_to_skip_decoy:
            native_residues.remove(residue)

    for counter,residue in enumerate(decoy.get_residues(),start=1):
        if counter in what_to_skip_native:
            decoy_residues.remove(residue)

    zipped_file = izip(native_residues,decoy_residues)


    if residues:
    # iterate through residues
        for(native_residue, decoy_residue) in zipped_file:
            res_native, res_decoy = native_residue.id[1], decoy_residue.id[1]
            chain_native, chain_decoy = native_residue.get_parent().id, decoy_residue.get_parent().id
            native_list = native_residue.get_list()
            decoy_list = decoy_residue.get_list()
            if (chain_native, res_native) in residues:
                try:
                    native_atoms_ca.append(native_residue['CA'])
                    decoy_atoms_ca.append(decoy_residue['CA'])
                except KeyError:
                    print "warning: residue", str(native_residue.get_id()[1]), "has no ca atom in either the native or the decoy structure.  either this residue is a hetatm or part of your backbone is missing. the residue will not be included in the rmsd calculations."
                try:
                    native_atoms_bb.append(native_residue['CA'])
                    native_atoms_bb.append(native_residue['C'])
                    native_atoms_bb.append(native_residue['O'])
                    native_atoms_bb.append(native_residue['N'])
                    decoy_atoms_bb.append(decoy_residue['CA'])
                    decoy_atoms_bb.append(decoy_residue['C'])
                    decoy_atoms_bb.append(decoy_residue['O'])
                    decoy_atoms_bb.append(decoy_residue['N'])
                except KeyError:
                    print "warning: residue", str(native_residue.get_id()[1]), "has no backbone atoms"

                zipped_atoms = izip(native_list, decoy_list)
                for (native_atom, decoy_atom) in zipped_atoms:
                    if native_atom.id != decoy_atom.id:
                        if debug:
                            print "atom {0} from chain {1}, residue {2} and atom {3} from chain {4}, residue {5} don't match up, skipping".format(
                                native_atom.id, chain_native, res_native, decoy_atom.id, chain_decoy, res_decoy)
                        continue
                    native_atoms_all.append(native_atom)
                    decoy_atoms_all.append(decoy_atom)
    # no residues specified
    else:
        for(native_residue, decoy_residue) in zipped_file:
            res_native, res_decoy = str(native_residue.id[1]), str(decoy_residue.id[1])
            chain_native, chain_decoy = native_residue.get_parent().id, decoy_residue.get_parent().id
            native_list = native_residue.get_list()
            decoy_list = decoy_residue.get_list()
            try:
                native_atoms_ca.append(native_residue['CA'])
                decoy_atoms_ca.append(decoy_residue['CA'])
            except KeyError:
                print "warning: residue", str(native_residue.get_id()[1]), "has no ca atom in either the native or the decoy structure.  either this residue is a hetatm or part of your backbone is missing. the residue will not be included in the rmsd calculations."
            try:
                native_atoms_bb.append(native_residue['CA'])
                native_atoms_bb.append(native_residue['C'])
                native_atoms_bb.append(native_residue['O'])
                native_atoms_bb.append(native_residue['N'])
                decoy_atoms_bb.append(decoy_residue['CA'])
                decoy_atoms_bb.append(decoy_residue['C'])
                decoy_atoms_bb.append(decoy_residue['O'])
                decoy_atoms_bb.append(decoy_residue['N'])
            except KeyError:
                print "warning: residue", str(native_residue.get_id()[1]), "has no backbone atoms"

            zipped_atoms = izip(native_list, decoy_list)
            for (native_atom, decoy_atom) in zipped_atoms:
                if native_atom.id != decoy_atom.id:
                    if debug:
                        print "atom {0} from chain {1}, residue {2} and atom {3} from chain {4}, residue {5} don't match up, skipping".format(
                            native_atom.id, chain_native, res_native, decoy_atom.id, chain_decoy, res_decoy)
                    continue
                native_atoms_all.append(native_atom)
                decoy_atoms_all.append(decoy_atom)

    superpose_ca = Superimposer()
    superpose_bb = Superimposer()
    superpose_all = Superimposer()

    superpose_ca.set_atoms(native_atoms_ca, decoy_atoms_ca)
    superpose_ca.apply(decoy.get_atoms())
    ca_rms = superpose_ca.rms
    superpose_bb.set_atoms(native_atoms_bb, decoy_atoms_bb)
    superpose_bb.apply(decoy.get_atoms())
    bb_rms = superpose_bb.rms
    superpose_all.set_atoms(native_atoms_all, decoy_atoms_all)
    superpose_all.apply(decoy.get_atoms())
    all_rms = superpose_all.rms

    return {"ca_rmsd": ca_rms, "bb_rmsd": bb_rms, "all_rmsd": all_rms}


def calculate_all_superpositions_clustal_constraints(native, decoy, residues, debug, native_list, decoy_list):
    """This basically uses mammoth to calculate which atoms should be thrown out of the calculation. See Mammoth LIB. It will return two lists\
    that you can iterate through at the same time and align each of those atoms. This will calculate all three RMSDS and can support a residue file from\
    the residue file class to just give RMSDS of certain segments of the PDB.
    The resfiles class first which returns a resfile list:\n
    1 G ALLAA\n
    3 H PIKAA G\n
    would return a list of tuples like :\n
    [(1,G),(3,H)]"""
    native_atoms = {}
    decoy_atoms = {}
    native_atoms_ca = []
    decoy_atoms_ca = []
    native_atoms_bb = []
    decoy_atoms_bb = []
    native_atoms_all = []
    decoy_atoms_all = []
    # the native and decoy list from mammoth class
    zipped_file = izip(native_list, decoy_list)

    # gather both decoy and native atoms
    for residue in native.get_residues():
        native_atoms[(residue.get_parent().id, int(residue.id[1]))] = {
            'atom_list': residue.get_list(), "name": residue.resname}

    for residue in decoy.get_residues():
        decoy_atoms[(residue.get_parent().id, int(residue.id[1]))] = {'atom_list': residue.get_list(), "name": residue.resname}

    if residues:
    # iterate through residues
        for(native_residue, decoy_residue) in zipped_file:
            if decoy_residue not in residues:
                continue
            if debug:
                print native_atoms[native_residue]['atom_list'][1].get_parent().resname, decoy_atoms[decoy_residue]['atom_list'][1].get_parent().resname
            for atom in native_atoms[native_residue]['atom_list']:
                if atom.id == 'CA':
                    native_atoms_ca.append(atom)
                    native_atoms_bb.append(atom)
                if atom.id == 'C' or atom.id == 'O' or atom.id == 'N':
                    native_atoms_bb.append(atom)
            # iterate through decoy atoms
            for atom in decoy_atoms[decoy_residue]['atom_list']:
                if atom.id == 'CA':
                    decoy_atoms_ca.append(atom)
                    decoy_atoms_bb.append(atom)
                if atom.id == 'C' or atom.id == 'O' or atom.id == 'N':
                    decoy_atoms_bb.append(atom)

            # one final check
            for native_atom, decoy_atom in zip(native_atoms[native_residue]['atom_list'], decoy_atoms[decoy_residue]['atom_list']):
                if native_atom.get_parent().resname != decoy_atom.get_parent().resname:
                    if native_atom.id == "CA" or native_atom.id == "O" or native_atom.id == "N" or native_atom.id == "C":
                        native_atoms_all.append(native_atom)
                    if decoy_atom.id == "CA" or decoy_atom.id == "O" or decoy_atom.id == "N" or decoy_atom.id == "C":
                        decoy_atoms_all.append(decoy_atom)
                    continue
                elif native_atom.id == decoy_atom.id:
                    native_atoms_all.append(native_atom)
                    decoy_atoms_all.append(decoy_atom)
    # no resfile specified
    else:
        for(native_residue, decoy_residue) in zipped_file:
            if debug:
                # print native_residue, decoy_residue
                print native_residue, native_atoms[native_residue]['name'], decoy_residue, decoy_atoms[decoy_residue]['name']
               # print native_atoms[native_residue][0].get_parent().resname, decoy_atoms[decoy_residue][0].get_parent().resname
            for atom in native_atoms[native_residue]['atom_list']:
                if atom.id == 'CA':
                    native_atoms_ca.append(atom)
                    native_atoms_bb.append(atom)
                if atom.id == 'C' or atom.id == 'O' or atom.id == 'N':
                    native_atoms_bb.append(atom)
            # iterate through decoy atoms
            for atom in decoy_atoms[decoy_residue]['atom_list']:
                if atom.id == 'CA':
                    decoy_atoms_ca.append(atom)
                    decoy_atoms_bb.append(atom)
                if atom.id == 'C' or atom.id == 'O' or atom.id == 'N':
                    decoy_atoms_bb.append(atom)

            # one final check
            for native_atom, decoy_atom in zip(native_atoms[native_residue]['atom_list'], decoy_atoms[decoy_residue]['atom_list']):
                if native_atom.get_parent().resname != decoy_atom.get_parent().resname:
                    if native_atom.id == "CA" or native_atom.id == "O" or native_atom.id == "N" or native_atom.id == "C":
                        native_atoms_all.append(native_atom)
                    if decoy_atom.id == "CA" or decoy_atom.id == "O" or decoy_atom.id == "N" or decoy_atom.id == "C":
                        decoy_atoms_all.append(decoy_atom)
                    continue
                elif native_atom.id == decoy_atom.id:
                    native_atoms_all.append(native_atom)
                    decoy_atoms_all.append(decoy_atom)

    superpose_ca = Superimposer()
    superpose_bb = Superimposer()
    superpose_all = Superimposer()
    superpose_ca.set_atoms(native_atoms_ca, decoy_atoms_ca)
    superpose_ca.apply(decoy.get_atoms())
    ca_rms = superpose_ca.rms

    superpose_bb.set_atoms(native_atoms_bb, decoy_atoms_bb)
    superpose_bb.apply(decoy.get_atoms())
    bb_rms = superpose_bb.rms

    superpose_all.set_atoms(native_atoms_all, decoy_atoms_all)
    superpose_all.apply(decoy.get_atoms())
    all_rms = superpose_all.rms

    return {"ca_rmsd": ca_rms, "bb_rmsd": bb_rms, "all_rmsd": all_rms}


def find_gaps(pdb, sequence, chain):
    """return a sequence with gaps in the pdb represented by - symbols"""
    # we can't trust the seqres record, it might not even exist, so get the sequence by looping through all the residues
    pdb_sequence = ""

    for residue in pdb.get_residues():
        if residue.get_full_id()[2] != chain and chain != 'ALL':  # skip residues that aren't in the chain
            continue
        residue_name = amino_acids.longer_names[residue.get_resname()]
        type(residue_name)
        pdb_sequence += residue_name
    if len(pdb_sequence) == 0:
        print "There is a problem with the command line options. Check your commandline for consitency with inputs."
    # now we align the two sequences
    alignment = pairwise2.align.globalmx(pdb_sequence, sequence, 2, -1)
    # print alignment
    return alignment[0][1]
