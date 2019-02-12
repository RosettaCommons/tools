#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   compare_info_files.py
## @brief  Script to compare difference between antibody database versions
## @author Jeliazko Jeliazkov

## @details compare two versions of the antibody database info files
## Loose note:
## In the old (non-auto) database, pdbs were named as "pdb123e_chothia.pdb"
## In the new, the naming convention is just "123e.pdb"
## Also, columns may differ

import collections
import glob
import os
import sys
from create_antibody_db import file_to_dict

def read_info_file(info_dir, regex_str):
    # helper function
    # store current dir for navigation
    cwd = os.getcwd()
    # go into dir
    os.chdir(os.path.abspath(info_dir))
    # e.g. we use either antibody.info or antibody_info
    # don't do this in the future... have user specify files
    possible_files = glob.glob(regex_str) # should only be one match

    if len(possible_files) > 1:
        print "Warning! Multiple possible files found!"
        print "Try a better regex than: " + regex_str
        print possible_files
        sys.exit()

    result = file_to_dict(possible_files[0])
    os.chdir(cwd)
    return result

def compare_ocds(overlap, old_ocds, new_ocds):
    # compare ocd values
    # note overlap for old file must be converted from four letter pdb
    # to pdbXXXX_chothia.pdb
    of = open("comparison/ocds.csv", "w")
    of.write("pdb1, pdb2, old_ocd, new_ocd\n")

    # comparison is across all pairs (yikes)
    miss_counts = 0

    for pdb1 in overlap:
        for pdb2 in overlap:
            if pdb1 != pdb2:
                pdb1_old = "pdb{}_chothia.pdb".format(pdb1)
                pdb2_old = "pdb{}_chothia.pdb".format(pdb2)

                OCDo = old_ocds[pdb1_old][pdb2_old]
                OCDn = new_ocds[pdb1][pdb2]

                if OCDo != OCDn:
                    of.write("{}, {}, {}, {}\n".format(pdb1, pdb2, OCDo, OCDn))
                    miss_counts += 1

    of.close()
    print "For OCDs, there were {} mismatches.".format(miss_counts)
    return

def compare_bfactors(cdr, overlap, old_bs, new_bs):
    # for overlapping pdbs, compare bfactors
    # new bfactors are reported as floats
    # old are T if > 50, else F
    of = open("comparison/{}_bfactors.csv".format(cdr), "w")
    of.write("cdr, pdb, old, new\n")

    miss_counts = 0

    for pdb in overlap:
        Bo = (old_bs[pdb][cdr] == "True")
        Bn = new_bs[pdb][cdr]

        if Bn != "none": Bn = (float(Bn) > 50.0)

        if Bo != Bn:
            of.write("{}, {}, {}, {}\n".format(cdr, pdb, Bo, Bn))
            miss_counts += 1

    of.close()
    print "For {} B factors, there were {} mismatches.".format(cdr, miss_counts)
    return

def compare_cdr(cdr, overlap, old_abi, new_abi):
    # for overlapping pdbs, spot check h1, h2, h3 ...
    of = open("comparison/{}_sequences.csv".format(cdr), "w")
    of.write("pdb, cdr, old, new\n")

    miss_counts = 0

    for pdb in overlap:
        if old_abi[pdb][cdr] != new_abi[pdb][cdr]:
            of.write("{}, {}, {}, {}\n".format(pdb, cdr, old_abi[pdb][cdr], new_abi[pdb][cdr]))
            miss_counts += 1

    of.close()
    print "For {} sequences, there were {} mismatches.".format(cdr, miss_counts)
    return

def compare_cdr_lens(cdr, old_cdri, new_cdri):
    # for all pdbs, spot check lengths
    of = open("comparison/{}_lengths.csv".format(cdr), "w")
    of.write("cdr, length, old, new\n")

    # hard coded for old manual file
    old_lens = [old_cdri[pdb]["{}_length".format(cdr.upper())] for pdb in old_cdri.keys()]
    old_lens = [int(x) for x in old_lens if x != "-"]
    old_counts = collections.Counter(old_lens)

    new_lens = [int(new_cdri[pdb]["{}_len".format(cdr)]) for pdb in new_cdri.keys()]
    new_counts = collections.Counter(new_lens)

    # write
    all_keys = old_counts.keys() + new_counts.keys()
    for key in all_keys:
        try:
            ov = old_counts[key]
        except KeyError:
            ov = 0
        try:
            nv = new_counts[key]
        except KeyError:
            nv = 0
        of.write("{}, {}, {}, {}\n".format(cdr, key, ov, nv))
    of.close()

    return

def compare(info_dir_new, info_dir_old):
    """
    Read files antibody, cdr, frh, frl, and frlh from directory.
    Compare columns pdb, resolution, h1, h2, h3, l1, l2, l3, frh, frl, light_heavy.
    Comparison does not have 100% coverage.
    """
    # compare antibody info files -- probably the most relevant
    old_abi = read_info_file(info_dir_old, "antibody*info")
    new_abi = read_info_file(info_dir_new, "antibody*info")

    # compare antibody info files -- probably the most relevant
    old_cdri = read_info_file(info_dir_old, "cdr*info")
    new_cdri = read_info_file(info_dir_new, "cdr*info")

    # compare bfactors
    old_bfactors = read_info_file("/Users/lqtza/Rosetta/main/database/protocol_data/antibody/", "list*bfactor50")
    new_bfactors = read_info_file(info_dir_new, "bfactor*info")

    # comapre OCDs
    old_ocds = read_info_file("/Users/lqtza/Rosetta/main/database/protocol_data/antibody/", "comparisons*txt")
    new_ocds = read_info_file(info_dir_new, "orientations*txt") # move to info file

    # print length of keys (pdbs) and overlap
    old_pdbs = old_abi.keys()
    new_pdbs = new_abi.keys()
    overlap_pdbs = list(set(old_pdbs) & set(new_pdbs))
    print "old, new, overlap"
    print "{}, {}, {}".format(len(old_pdbs), len(new_pdbs), len(overlap_pdbs))

    # comparisons write to comparison directory
    if not os.path.isdir("comparison"): os.mkdir("comparison")

    # compare sequences
    compare_cdr("h1", overlap_pdbs, old_abi, new_abi)
    compare_cdr("h2", overlap_pdbs, old_abi, new_abi)
    compare_cdr("h3", overlap_pdbs, old_abi, new_abi)
    compare_cdr("l1", overlap_pdbs, old_abi, new_abi)
    compare_cdr("l2", overlap_pdbs, old_abi, new_abi)
    compare_cdr("l3", overlap_pdbs, old_abi, new_abi)

    # compare length distributions (no need for overlap
    compare_cdr_lens("h1", old_cdri, new_cdri)
    compare_cdr_lens("h2", old_cdri, new_cdri)
    compare_cdr_lens("h3", old_cdri, new_cdri)
    compare_cdr_lens("l1", old_cdri, new_cdri)
    compare_cdr_lens("l2", old_cdri, new_cdri)
    compare_cdr_lens("l3", old_cdri, new_cdri)

    # compare b factors -- repeat overlap calc ???
    # print length of keys (pdbs) and overlap
    old_pdbs = old_bfactors.keys()
    new_pdbs = new_bfactors.keys()
    overlap_pdbs = list(set(old_pdbs) & set(new_pdbs))
    print "for bfactor file!"
    print "old, new, overlap"
    print "{}, {}, {}".format(len(old_pdbs), len(new_pdbs), len(overlap_pdbs))
    compare_bfactors("h1", overlap_pdbs, old_bfactors, new_bfactors)
    compare_bfactors("h2", overlap_pdbs, old_bfactors, new_bfactors)
    compare_bfactors("h3", overlap_pdbs, old_bfactors, new_bfactors)
    compare_bfactors("l1", overlap_pdbs, old_bfactors, new_bfactors)
    compare_bfactors("l2", overlap_pdbs, old_bfactors, new_bfactors)
    compare_bfactors("l3", overlap_pdbs, old_bfactors, new_bfactors)


    # compare OCDs -- repeat overlap calc ???
    # convert pdbXXXX_chothia.pdb into XXXX
    old_pdbs = [x[3:7] for x in old_ocds.keys()]
    new_pdbs = new_ocds.keys()
    overlap_pdbs = list(set(old_pdbs) & set(new_pdbs))
    print "for OCD file!"
    print "old, new, overlap"
    print "{}, {}, {}".format(len(old_pdbs), len(new_pdbs), len(overlap_pdbs))
    compare_ocds(overlap_pdbs, old_ocds, new_ocds)

    return

if __name__ == "__main__":
    compare("/Users/lqtza/Rosetta/tools/antibody-update/info/", "/Users/lqtza/Rosetta/tools/antibody/info/")
    #compare("/Users/jjeliazkov/Rosetta/tools/antibody-update/info/", "/Users/jjeliazkov/Rosetta/tools/antibody/info/")
