#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   create_antibody_db.py
## @brief  Script to create database for RosettaAntibody
## @author Jeliazko Jeliazkov

## @details download non-redundant Chothia Abs from SAbDab
## Abs are downloaded by html query (is there a better practice)?
## Abs are Chothia-numbered, though we use Kabat to define CDRs.
## After download, trim Abs to Fv and extract FR and CDR sequences.
## For the purpose of trimming, truncate the heavy @112 and light @109
## Some directories (antibody_database, info, etc...) are hard coded.

import urllib2
import math
import re
import os
import bz2
import sys
from collections import defaultdict
from subprocess import call

def get_cdr_ranges():
    """
    Central defintion for CDR ranges.
    """

    cdr_ranges = {  'h1': [26, 35], # res 31 has icode
                    'h2': [50, 65], # res 52 has icode
                    'h3': [95, 102], # res 95 has icode
                    'l1': [24, 34], # res 30 has icode
                    'l2': [50, 56],
                    'l3': [89, 97] } # res 95 has icode
    return cdr_ranges

def get_fr_ranges():
    """
    Central definition for FR ranges, inclusive as above.
    """
    fr_ranges = { 'frh': list(range(10, 26)) + list(range(36, 40)) + list(range(46, 50)) + list(range(66, 95)) + list(range(103, 110)),
                  'frl': list(range(10, 24)) + list(range(35, 40)) + list(range(46, 50)) + list(range(57, 67)) + list(range(71, 89)) + list(range(98, 105)) }

    return fr_ranges

def check_pdb_line_for_atoms(line, atoms):
    present = False
    for atom in atoms:
        if line[12:16].strip() == atom:
            present = True
    return present

def get_avg_bfactor_from_pdb(pdb_text, chain, res_range):
    """ Given pdb text, chain, and resnum range, <B>. """
    B = 0.0
    natoms = 0
    for line in pdb_text.split("\n"):
        #if line.startswith("ATOM") and line[12:16].strip() == "CA": # if only CA matters
        if line.startswith("ATOM"): # all atoms
        #if line.startswith("ATOM") and check_pdb_line_for_atoms(line, ["N", "CA", "C", "O"]):
            cres = int(line[22:26])
            if line[21] == chain and cres >= res_range[0] and cres <= res_range[1]:
                B += float(line[60:66])
                natoms += 1

    # sometimes loops are not present, return "none" for these cases?
    if natoms == 0: return "none"

    return B/natoms

def get_bfactors():
    """
    Compute average Bs for all CDRs.
    Write to file.
    """

    headers = ['h1', 'h2', 'h3', 'l1', 'l2', 'l3']
    cdr_ranges = get_cdr_ranges() # same headers as above
    bfactors = {}

    # iterate over PDBs in antibody_database and truncate accordingly
    unique_pdbs = set([x[:4] for x in os.listdir("antibody_database") if x.endswith(".pdb") or x.endswith(".pdb.bz2")])

    for pdb in unique_pdbs:
        pdb_text = ""
        try:
            with open("antibody_database/" + pdb + ".pdb", "r") as f:
                pdb_text = f.read() # want string not list
        except IOError:
            sys.exit("Failed to open {} in antibody_database/ !".format(pdb))

        # should have pdb_text now
        if len(pdb_text) == 0: sys.exit("Nothing parsed for PDB {} !".format(pdb))

        bfactors[pdb] = {"pdb":pdb} # initialize results dict
        for col in headers:
            bfactors[pdb][col] = get_avg_bfactor_from_pdb(pdb_text, col[0].upper(), cdr_ranges[col])

    dict_to_file("info/bfactors.info", ["pdb"] + headers, bfactors)

    return

def calc_orientation_diff(o1, o2):
    """
    Calculated as a z-score. Mean/Stdev are from Nick's paper:
    https://academic.oup.com/peds/article/29/10/409/2462315
    See Figure 3
    """
    # order is distance, Hopen, Lopen, packing_angle
    # mimics output of vl_vh_orientation_coords function
    means = [14.6, 97.2, 99.4, -52.3]
    sigmas = [0.32, 2.55, 1.93, 3.83]
    res = 0.0
    for i in [1,2,3,4]: # Rosetta indexing :O why ?
        res += ((o1[i] - o2[i])/sigmas[i-1])**2
    return res

def compare_orientations():
    """
    Necessary for multitemplate grafting.
    Use Nick's pilot app (or PyRosetta) to compute LHOCs and write to file.
    """
    # last to do!! SLOW
    # could speed up by moving Nick's calculation in here and assuming Chothia numbering
    # place results in: protocol_data/antibody/comparisons.txt
    # or update to place data in info

    # loop over all *paired* PDBs, load and calculate orientations, then calculate OCDs (maybe dump halfway?)
    frlh_info = file_to_dict("info/frlh.info")
    unique_pdbs = frlh_info.keys()
    # d[pdb] = [distance, heavy_opening_angle, light_opening_angle, packing_angle]
    orientations = {}
    # for fun sum all the values and calculate averages + stdevs
    orientation_sums = [0.0, 0.0, 0.0, 0.0]

    for pdb in unique_pdbs:

        print("Reading for OCD calculations: " + pdb + "...")

        pose = pyrosetta.pose_from_file("antibody_database/" + pdb + ".pdb")
        abi = antibody.AntibodyInfo(pose)

        orientations[pdb] = antibody.vl_vh_orientation_coords(pose, abi)
        # check that orientations aren't all zero and if so, warn!
        if sum(orientations[pdb]) == 0:
            print("Warning read PDB, but couldn't calculate OCDs!")
            print("Removing {} from antibody_database/".format(pdb))
            os.remove("antibody_database/" + pdb + ".pdb")
            del orientations[pdb]
            sys.exit("Please re-run all steps *after* truncate.")
        # list iteration to sum
        orientation_sums = [sum(x) for x in zip(orientation_sums, orientations[pdb])]

    # write all the values for the individual pdbs
    with open("info/angles.info", "w") as of:
        of.write("pdb, distance, heavy_opening_angle, light_opening_angle, packing_angle\n")
        for pdb in unique_pdbs: # break writing paradigm -- sorry!
            of.write("{}, {:.3f}, {:.3f}, {:.3f}, {:.3f}\n".format(pdb, *(orientations[pdb])))

    # for fun, print averages over database
    print("Average relative orientation values (+/-) a standard deviation: ")
    for i, meas in enumerate(["distance", "heavy_opening_angle", "light_opening_angle", "packing_angle"]):
        mm = orientation_sums[i]/len(unique_pdbs)
        mstdev = 0.0
        for pdb in unique_pdbs:
            mstdev += (orientations[pdb][i+1] - mm)**2 # Rosetta indexing
        mstdev = math.sqrt(mstdev/len(unique_pdbs))
        print("\t {}: {:.2f} +/- {:.2f}".format(meas, mm, mstdev))
    print("Originally reported by Marze et al. in PEDS (2016).")

    # loop over all pairs
    pairs = {}
    for p1 in unique_pdbs:
        pairs[p1] = {"pdb":p1}
        for p2 in unique_pdbs:
            pairs[p1][p2] = calc_orientation_diff(orientations[p1], orientations[p2])

    dict_to_file("info/orientations.txt", ["pdb"] + unique_pdbs, pairs)
    return

def tuple_list_to_fasta(tl, fn):
    """
    Take a list of pdb, sequence pairs and output a single fasta file.
    Input tuple must be (pdb, seq)
    """
    with open(fn, "w") as f:
        for pair in tl:
            f.write(">" + pair[0] + "\n")
            f.write(pair[1] + "\n")
    return

def create_blast_db():
    """
    Using info files, extract sequences per length and setup blast databases.
    We need one database per CDR and length pair, as well as:
    FRH, FRL, light_heavy (rename to orientation) databases.
    """

    cdr_info = file_to_dict("info/cdr.info")

    # convert to dict of dicts set up as d[CDR][Length] = [(pdb, "fasta string"),...]
    converted_dict = { "h1":defaultdict(list),
                        "h2":defaultdict(list),
                        "h3":defaultdict(list),
                        "l1":defaultdict(list),
                        "l2":defaultdict(list),
                        "l3":defaultdict(list)}

    for pdb in cdr_info.keys():
        for cdr in converted_dict.keys():
            try:
                length = cdr_info[pdb][cdr+"_len"]
                converted_dict[cdr][length].append((pdb, cdr_info[pdb][cdr]))
            except KeyError: # length unmatched
                continue

    if not os.path.isdir("blast_database"): os.mkdir("blast_database")

    for cdr in converted_dict.keys():
        for length in converted_dict[cdr].keys():
            tuple_list_to_fasta(converted_dict[cdr][length], "blast_database/database.{}.{}".format(cdr.upper(),length))

            call(["makeblastdb",
                    "-in", "blast_database/database.{}.{}".format(cdr.upper(),length),
                    "-dbtype", "prot",
                    "-title", "database.{}.{}".format(cdr.upper(),length),
                    "-out", "blast_database/database.{}.{}".format(cdr.upper(),length)])

    # next code for the FRH/FRL/orientation databases
    frh_info = file_to_dict("info/frh.info") # same as heavy?
    frl_info = file_to_dict("info/frl.info") # same as light?
    frlh_info = file_to_dict("info/frlh.info") # referred to as light_heavy previously

    # no need for converted dict, only one structural region per file
    # so extract tuples of pdb, sequence from the correct columns
    frh_tuples = []
    frl_tuples = []
    frlh_tuples = []

    for pdb in frh_info.keys(): frh_tuples.append((pdb, frh_info[pdb]["frh"]))
    for pdb in frl_info.keys(): frl_tuples.append((pdb, frl_info[pdb]["frl"]))
    for pdb in frlh_info.keys(): frlh_tuples.append((pdb, frlh_info[pdb]["light_heavy"]))

    # now write fastas, then make dbs
    tuple_list_to_fasta(frh_tuples, "blast_database/database.FRH")
    tuple_list_to_fasta(frl_tuples, "blast_database/database.FRL")
    tuple_list_to_fasta(frlh_tuples, "blast_database/database.light_heavy")

    call(["makeblastdb",
            "-in", "blast_database/database.FRH",
            "-dbtype", "prot",
            "-title", "database.FRH",
            "-out", "blast_database/database.FRH"])

    call(["makeblastdb",
            "-in", "blast_database/database.FRL",
            "-dbtype", "prot",
            "-title", "database.FRL",
            "-out", "blast_database/database.FRL"])

    call(["makeblastdb",
            "-in", "blast_database/database.light_heavy",
            "-dbtype", "prot",
            "-title", "database.light_heavy",
            "-out", "blast_database/database.light_heavy"])

    return

def file_to_dict(fn):
    """
    Very specific function reads output of func below.
    Relies on a header "# pdb col1 col2 ...".
    Parses into dict[pdb] = [col1: val1, col2: val2]
    """

    d = {}

    with open(fn, "r") as f:

        header = f.readline()
        if not header.startswith("#"): sys.exit("Invalid file: {}!".format(fn))
        header = header.strip().split()[2:]

        for line in f.readlines():

            sl = line.split()
            pdb = sl[0]
            d[pdb] = {}

            for col,val in zip(header,sl[1:]):
                d[pdb][col] = val

    return d

def dict_to_file(fn, header, info_dict):
    """
    Very specific function to write dict of dicts to file.
    Keys for first dict are pdbs.
    Keys for second dict are column names.
    """
    with open(fn, "w") as f:
        # write header
        f.write("# " + ' '.join(header) + "\n")
        for pdb in info_dict.keys():
            line = ""
            for column in header:
                if info_dict[pdb][column] == "":
                    line += " none"
                else:
                    line += " {}".format(info_dict[pdb][column])
            f.write(line.strip() + "\n")
    return

def extract_cdr_sequence(pose, chain, start, stop):
    """
    Go from start to stop residue. Check for loop integrity along the way.
    If the loop is intact, then the sequence must be good.
    """
    sequence = ""
    pose_start = pose.pdb_info().pdb2pose(chain, start)
    pose_stop = pose.pdb_info().pdb2pose(chain, stop)

    if pose_start == 0 or pose_stop == 0:
        return ""

    for i in range(pose_start, pose_stop+1):
        has_issues = loops.has_severe_pep_bond_geom_issues(pose, i)
        # geometry issues -- skip
        if has_issues[0]: return ""
        sequence += pose.sequence()[i-1]

    return sequence

def extract_fr_sequence(pose, chain, residue_list):
    """
    Extract the sequence for a certain region of a Chothia-numbered structure.
    Will break for regions shorter than expected
    Sequence will be length 0 if the region is missing residues or
    if the bond lengths/angles fail a geometry check.

    Careful with insertions codes around residue 82

    Otherwise, iterate through the list and convert.
    """
    sequence = ""
    rosetta_numbers = []
    # loop over certain regions and extract seq / check quality
    # first convert chothia numbering to rosetta numbering
    for i in residue_list:

        pose_num = pose.pdb_info().pdb2pose(chain, i)

        # missing residue, skip region
        if pose_num == 0: return ""

        rosetta_numbers.append(pose_num)
        # account for icodes by inserting
        if chain == "H" and i == 82:
            # add pose numbers until the next residue
            rosetta_numbers.append(pose_num + 1) #82A
            rosetta_numbers.append(pose_num + 2) #82B
            rosetta_numbers.append(pose_num + 3) #82C


    for i in rosetta_numbers:
        # residue cannot be terminal because the function below checks +1
        # structure terminated mid region in the case below
        if i == pose.size(): return ""

        has_issues = loops.has_severe_pep_bond_geom_issues(pose, i)
        # geometry issues, skip region
        if has_issues[0]: return ""

        # all good -- append loop sequence
        sequence += pose.sequence()[i-1]

    return sequence

def write_info_files():
    """
    Several quality filters need info files.
    Generate these from the structures + SAbDab summary file.
    """

    # loop over all PDBs in the database and match them to SAbDab records
    sabdab_dict = parse_sabdab_summary("info/sabdab_summary.tsv")

    # iterate over PDBs in antibody_database and truncate accordingly
    unique_pdbs = set([x[:4] for x in os.listdir("antibody_database") if x.endswith(".pdb") or x.endswith(".pdb.bz2")])

    # dicts of info to store
    infos = {'antibody': {},
                'cdr': {},
                'frh': {},
                'frl': {},
                'frlh': {}}

    # different info files store different info...
    # use SAbDab headers for lazy matching
    headers = {'antibody': ['pdb', 'resolution', 'date', 'light_ctype', 'h1', 'h2', 'h3', 'l1', 'l2', 'l3', 'frh', 'frl'],
                'cdr': ['pdb', 'resolution', 'date', 'light_ctype', 'h1', 'h2', 'h3', 'l1', 'l2', 'l3', 'h1_len', 'h2_len', 'h3_len', 'l1_len', 'l2_len', 'l3_len'],
                'frh': ['pdb', 'resolution', 'heavy_len', 'frh_len', 'frh'],
                'frl': ['pdb', 'resolution', 'light_len', 'frl_len', 'frl'],
                'frlh': ['pdb', 'light_heavy_len', 'light_heavy']}

    # residue ranges for cdrs
    cdr_ranges = get_cdr_ranges()
    # could move to function if used elsewhere in the code
    fr_ranges =  get_fr_ranges()

    # frh missing conserved residues
    missing_frh_cons = []
    # frl missing conserved resiudes
    missing_frl_cons = []

    for pdb in unique_pdbs:
        print("extracting info from " + pdb + "...")
        # store features under column names for each pdb (easy matching)
        features = {}
        features['pdb'] = pdb
        # easy features from summary file
        features['resolution'] = sabdab_dict[pdb]['resolution']
        features['date'] = sabdab_dict[pdb]['date']
        features['light_ctype'] = sabdab_dict[pdb]['light_ctype']
        # ok now more challenging features such as frh/frl/cdr sequences
        # these come from the structure
        try:
            with open("antibody_database/" + pdb + ".pdb", "r") as f:
                pdb_text = f.read() # want string not list
        except IOError:
            sys.exit("Failed to open {} in antibody_database/ !".format(pdb))

        # if we are loading in Rosetta then we should use pose2pdb to get the
        # sequence as this would help us establish missing residues
        # do some quality checks on geometry
        # we don't want to graft something that has a chain break
        # strong assumption here that we can load the file into Rosetta

        pose = pyrosetta.pose_from_file("antibody_database/" + pdb + ".pdb")

        features["h1"] = extract_cdr_sequence(pose, "H", *cdr_ranges["h1"])
        features["h2"] = extract_cdr_sequence(pose, "H", *cdr_ranges["h2"])
        features["h3"] = extract_cdr_sequence(pose, "H", *cdr_ranges["h3"])

        features["l1"] = extract_cdr_sequence(pose, "L", *cdr_ranges["l1"])
        features["l2"] = extract_cdr_sequence(pose, "L", *cdr_ranges["l2"])
        features["l3"] = extract_cdr_sequence(pose, "L", *cdr_ranges["l3"])

        # note fr ranges are full lists, not just start/stop as above
        features["frh"] = extract_fr_sequence(pose, "H", fr_ranges["frh"])
        features["frl"] = extract_fr_sequence(pose, "L", fr_ranges["frl"])

        # get light/heavy paired seq for orientation -- structural flaws ok here
        # hence we pull from the pdb_text directly rather than Rosetta
        features["light_heavy"] = get_sequence_from_pdb(pdb_text, "L", [5, 104])
        features["light_heavy"] += get_sequence_from_pdb(pdb_text, "H", [5, 109])


        # ranges from protocols/antibody/grafting/grafter.cc
        # compact representation
        conserved_frh_residues = list(range(10,26))
        conserved_frh_residues += list(range(36,50))
        conserved_frh_residues += [66]
        conserved_frh_residues += list(range(69,95))
        conserved_frh_residues += list(range(103,106))

        conserved_frl_residues = list(range(10,24))
        conserved_frl_residues += list(range(35,50))
        conserved_frl_residues += list(range(57,67))
        conserved_frl_residues += list(range(69,89))
        conserved_frl_residues += list(range(98,101))

        # conserved residues are necessary to superipmose FRH/L models
        # onto orientation model, so we need an extract check to exclude
        # frh/frl/orientation models with missing conserved residues
        if extract_fr_sequence(pose, "H", conserved_frh_residues) == "":
            # missing residue or bad geometry, exclude frh + orientation
            features["frh"] = ""
            features["frh_len"] = 0
            features["light_heavy"] = ""
            missing_frh_cons.append(pdb)

        if extract_fr_sequence(pose, "L", conserved_frl_residues) == "":
            # missing residue or bad geometry, exclude frh + orientation
            features["frl"] = ""
            features["frl_len"] = 0
            features["light_heavy"] = ""
            missing_frl_cons.append(pdb)

        # try to load the antibody info object
        # if we can't then we can't calculate OCD, so this Ab should not be
        # used as an orientation template
        try:
            abi = antibody.AntibodyInfo(pose)
        except: # if Ab cannot be constructed, something terrible is happening
            features["light_heavy"] = ""

        # set lengths last
        # get lengths for heavy
        features["h1_len"] = len(features["h1"])
        features["h2_len"] = len(features["h2"])
        features["h3_len"] = len(features["h3"])

        features["frh_len"] = len(features["frh"])

        features["l1_len"] = len(features["l1"])
        features["l2_len"] = len(features["l2"])
        features["l3_len"] = len(features["l3"])

        features["frl_len"] = len(features["frl"])

        features["light_heavy_len"] = len(features["light_heavy"])

        # not 100% sure about the custom length here...
        # also, we don't really blast against this so, let's just pull from text
        # in the future maybe remove this information?
        features["light_len"] = len(get_sequence_from_pdb(pdb_text, "L", [5, 104]))
        features["heavy_len"] = len(get_sequence_from_pdb(pdb_text, "H", [5, 109]))

        # store data, but how to handle "" ?
        # set up antibody_info
        td = {}
        for column in headers["antibody"]:
            td[column] = features[column]
        infos["antibody"][pdb] = td

        # set up cdr_info
        td = {}
        for column in headers["cdr"]:
            td[column] = features[column]
        infos["cdr"][pdb] = td

        # set up frh_info, only if there is an FRH
        if features["frh_len"] > 0:
            td = {}
            for column in headers["frh"]:
                td[column] = features[column]
            infos["frh"][pdb] = td

        # set up frl_info, only if there is an FRL
        if features["frl_len"] > 0:
            td = {}
            for column in headers["frl"]:
                td[column] = features[column]
            infos["frl"][pdb] = td

        # set up frlh_info, only if there are both chains
        # double check these lengths
        if features["light_heavy_len"] > 0:
            td = {}
            for column in headers["frlh"]:
                td[column] = features[column]
            infos["frlh"][pdb] = td

    # done looping over PDBs

    # report frh/frl issues
    print("{} pdbs had missing conserved_frh_residues!".format(len(missing_frh_cons)))
    print("{} pdbs had missing conserved_frl_residues!".format(len(missing_frl_cons)))

    # write different dicts to files
    for key in infos.keys():
        print("writing info/{}.info file".format(key))
        dict_to_file("info/{}.info".format(key), headers[key], infos[key])

    return

def parse_sabdab_summary(file_name):
    """
    SAbDab produces a rather unique summary file.
    This function reads that file into a dict with the key being the
    4-letter PDB code.
    """
    # dict[pdb] = {col1 : value, col2: value, ...}
    sabdab_dict = {}

    with open(file_name, "r") as f:
        # first line is the header, or all the keys in our sub-dict
        header = f.readline().strip().split("\t")
        # next lines are data
        for line in f.readlines():
            split_line = line.strip().split("\t")
            td = {} # temporary dict of key value pairs for one pdb
            for k,v in zip(header[1:], split_line[1:]):
                # pdb id is first, so we skip that for now
                td[k] = v
            # add temporary dict to sabdab dict at the pdb id
            sabdab_dict[split_line[0]] = td

    return sabdab_dict

def get_sequence_from_pdb(pdb_text, chain, res_range):
    """ Given pdb text, chain, and resnum range, return sequence. """
    seq = ""
    # dict to convert three residue code to one -- probably exists elsewhere...
    aa3to1 = {  'ALA':'A',
                'CYS':'C',
                'ASP':'D',
                'GLU':'E',
                'PHE':'F',
                'GLY':'G',
                'HIS':'H',
                'ILE':'I',
                'LYS':'K',
                'LEU':'L',
                'MET':'M',
                'ASN':'N',
                'PRO':'P',
                'GLN':'Q',
                'ARG':'R',
                'SER':'S',
                'THR':'T',
                'VAL':'V',
                'TRP':'W',
                'TYR':'Y'}
    for line in pdb_text.split("\n"):
        if line.startswith("ATOM") and line[12:16].strip() == "CA": # only CA matters
            cres = int(line[22:26])
            if line[21] == chain and cres >= res_range[0] and cres <= res_range[1]:
                seq += aa3to1[line[17:20]]

    return seq

def truncate_chain(pdb_text, chain, resnum, newchain):
    """
    Read PDB line by line and return all lines for a chain,
    with a resnum less than or equal to the input.
    This has to be permissive for insertion codes.
    This will return only a single truncated chain.
    This function can update chain to newchain.
    """
    trunc_text = ""
    for line in pdb_text.split("\n"):
        if line.startswith("ATOM") and line[21] == chain and int(line[22:26]) <= resnum:
            trunc_text += line[:21] + newchain + line[22:]
            trunc_text += "\n"
    return trunc_text

def truncate_antibody_pdbs():
    """
    We only use the Fv as a template, so this function loads each pdb
    and deletes excess chains/residues. We define the Fv, under the
    Chothia numbering scheme as H1-H112 and L1-L109.
    """
    # count warnings
    warn_pdbs = []

    # count delete files (without sabdab info)
    remove_pdbs = []

    # count deleted files (VH+VL same chain)
    same_chain = []

    # read SAbDab info for chain identities
    sabdab_dict = parse_sabdab_summary("info/sabdab_summary.tsv")

    # iterate over PDBs in antibody_database and truncate accordingly
    #unique_pdbs = set([x[:4] for x in os.listdir("antibody_database") if x.endswith(".pdb") or x.endswith(".pdb.bz2")])
    unique_pdbs = set([x[:4] for x in os.listdir("antibody_database") if x.endswith(".pdb")])

    for pdb in unique_pdbs:
        # try reading bzipped pdb, then regular pdb
        print("Truncating " + pdb + "...")
        pdb_text = ""
        try:
            with open("antibody_database/" + pdb + ".pdb", "r") as f:
                pdb_text = f.read() # want string not list
        except IOError:
            sys.exit("Failed to open {} in antibody_database/ !".format(pdb))

        # should have pdb_text now
        if len(pdb_text) == 0: sys.exit("Nothing parsed for PDB {} !".format(pdb))

        # test if pdb is in sabdab summary file, if not skip and delete PDB from db
        try:
            sabdab_dict[pdb]
        except KeyError:
            remove_pdbs.append(pdb)
            print(pdb + " not in sabdab summary file, removing ...")
            os.remove("antibody_database/" + pdb + ".pdb")
            continue

        hchain = sabdab_dict[pdb]["Hchain"]
        hchain_text = ""
        lchain = sabdab_dict[pdb]["Lchain"]
        lchain_text = ""

        # we do not currently have a good way of handling VH & VL on the same chain
        if hchain == lchain:
            same_chain.append(pdb)
            print(pdb + " has the VH+VL on a single chain, removing...")
            os.remove("antibody_database/" + pdb + ".pdb")
            continue

        if not hchain == "NA":
            hchain_text = truncate_chain(pdb_text, hchain, 112, "H")
            if len(hchain_text) == 0:
                # could not find heavy chain -- do not overwrite, but warn!
                warn_pdbs.append(pdb)
                print("Warning, could not find " + hchain + " chain for " + pdb + " !")
                print("It was not reported to be NA, so the file may have been altered!")
                continue

        if not lchain == "NA":
            lchain_text = truncate_chain(pdb_text, lchain, 109, "L")
            if len(lchain_text) == 0:
                # could not find heavy chain -- do not overwrite, but warn!
                warn_pdbs.append(pdb)
                print("Warning, could not find " + lchain + " chain for " + pdb + " !")
                print("It was not reported to be NA, so the file may have been altered!")
                continue

        # overwrite -- dangerous?
        with open("antibody_database/" + pdb+".pdb", "w") as f:
            f.write(hchain_text + lchain_text)

    for pdb in remove_pdbs:
        print("Deleted " + pdb + " from database because it is missing from summary file")

    for pdb in same_chain:
        print("Deleted " + pdb + " from database because it has VH+VL on the same chain.")
    print("Deleted {} total of same chain VH+VLs.".format(len(same_chain)))

    if len(warn_pdbs) > 0:
        print("Finished truncating, with {} warnings.".format(len(warn_pdbs)))
        #sys.exit("Exiting prematurely due to warnings.")

    return

def download_antibody_pdbs():
    """
    Function queries SAbDab for the latest list of non-redundant,
    chothia-numbered, sub 3-Angstrom antibody structures.
    """
    # base url for SAbDab query
    url_base = 'http://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/'
    url_query_base = url_base + 'DBrowser.php?'

    # configure options for SAbDab query
    options = { 'submitted' : 'true',
                'nonredundant' : 'true',
                'ab_ident' : '99%25',# 25 need for html
                'only_complete' : 'false', # might need capital F
                'incomple_nr' : '*',
                'agtype' : '*',
                'ag_ident' : '99%25',
                'aglengthnr' : '50',
                'resolution' : '3.0',
                'nonredundant_selected' : 'true' }

    # alternative options (all sub 3.0A structures)
    options2 = { 'submitted' : 'true',
                 'advanced' : 'true',
                 'abtype' : '*',
                 'method' : '*',
                 'organism' : '*',
                 'resolution' : '3.0',
                 'rfactor': '*',
                 'incomplex' : '*',
                 'ag_type' : '*',
                 'ag_length' : '50',
                 'lctype' : '*',
                 'hasconstant' : '*',
                 'hasaffinity' : '*' }

    # generate url from query
    #sabdab_url =  url_query_base + '&'.join(['{0}={1}'.format(x,y) for x,y in options.items()])
    sabdab_url =  url_query_base + '&'.join(['{0}={1}'.format(x,y) for x,y in options2.items()])

    # open page
    webpage = urllib2.urlopen(urllib2.Request(sabdab_url))

    # read html and extract summary file, and download chothia PDBs
    webpage_html = webpage.read() # PDBs in table, summary at end

    # use regex to extract PDBs
    chothia_pdb_re = re.compile('data/entries/\w\w\w\w/structure/chothia/\w\w\w\w.pdb')
    summary_re = re.compile('a href=.* download=\"summary.csv\"')
    pdb_urls = chothia_pdb_re.findall(webpage_html, re.MULTILINE)
    summary_url = summary_re.findall(webpage_html, re.MULTILINE) # should only find two
    summary_url = summary_url[0].split(" download")[0][8:-1] # drop unnecessary chars

    # make info and antibody_database directories if not present
    if not os.path.isdir("info"): os.mkdir("info")
    if not os.path.isdir("antibody_database"): os.mkdir("antibody_database")

    # download summary file, always overwrite
    with open("info/sabdab_summary.tsv", "w") as f:
        u = urllib2.urlopen(urllib2.Request(url_base + summary_url))
        f.write(u.read())

    # download PDBs, check for *.pdb or *.pdb.bz before writing
    print("Found {} pdbs, beginning download.".format(len(pdb_urls)))
    counter = 0
    for pdb in pdb_urls:
        if pdb[-8:-4] in ["1qok", "1dzb", "6b0w"]: # bad geometry -- skip
            continue
        if pdb[-8:-4] in ["6db7", "6iut"]: # AntibodyInfo construction issue
            continue
        counter += 1
        fpath = "antibody_database/{}.pdb".format(pdb[-8:-4])
        #fpath_bz = "antibody_database/{}.pdb.bz2".format(pdb[-8:-4])
        #if not (os.path.isfile(fpath) or os.path.isfile(fpath_bz)):
        if not os.path.isfile(fpath):
            with open(fpath, "w") as f:
                print("Downloading {}... {}/{}".format(pdb[-8:-4], counter, len(pdb_urls)))
                # won't timeout if no internet, so ...
                u = urllib2.urlopen(urllib2.Request(url_base + pdb))
                # write pdb ?
                f.write(u.read())
                # write compressed bz2
                #f.write(bz2.compress(u.read()))

    return

def create_antibody_db():
    """
    Run all function require to setup the database.
    """
    #download_antibody_pdbs()
    #truncate_antibody_pdbs()
    write_info_files()
    #create_blast_db()
    #compare_orientations()
    #get_bfactors()
    return

if __name__ == "__main__":
    # we'll eventually need Rosetta, so let's import it now (once)
    import pyrosetta
    from pyrosetta.rosetta.protocols import antibody
    from pyrosetta.rosetta.protocols import loops
    pyrosetta.init("-check_cdr_chainbreaks false")
    # check for execution in correct dir
    # or risk creation of dirs/files elsewhere
    if not os.getcwd().partition("tools/")[2] == "antibody-update":
        sys.exit("script needs to be run in tools/antibody-update!")
    create_antibody_db()
