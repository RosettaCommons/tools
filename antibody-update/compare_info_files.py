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
    possible_files = glob.glob(regex_str) # should only be one match

    if len(possible_files) > 1: 
        print "Warning! Multiple possible files found!"
        print "Try a better regex than: " + regex_str
        print possible_files
        sys.exit()

    result = file_to_dict(possible_files[0])
    os.chdir(cwd)
    return result

def compare_cdr(cdr, overlap, old_abi, new_abi):
    # for overlapping pdbs, spot check h1, h2, h3 ...
    print "COMPARING {}s!".format(cdr)
    print "\t{} mismatch: pdb, old, new".format(cdr)
    for pdb in overlap:
        if old_abi[pdb][cdr] != new_abi[pdb][cdr]:
            print "\t{} mismatch: {}, {}, {}".format(cdr, pdb, old_abi[pdb][cdr], new_abi[pdb][cdr])
    return

def compare_cdr_lens(cdr, old_cdri, new_cdri):
    # for all pdbs, spot check lengths
    print "COMPARING {}s!".format(cdr)

    # hard coded for old manual file
    old_lens = [old_cdri[pdb]["{}_length".format(cdr.upper())] for pdb in old_cdri.keys()]
    old_lens = [int(x) for x in old_lens if x != "-"]
    old_counts = collections.Counter(old_lens)

    new_lens = [int(new_cdri[pdb]["{}_len".format(cdr)]) for pdb in new_cdri.keys()]
    new_counts = collections.Counter(new_lens)

    print "OLD"
    print old_counts.keys()
    print old_counts.values()
    print ""

    print "New"
    print new_counts.keys()
    print new_counts.values()
    print ""

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

    # print length of keys (pdbs) and overlap
    old_pdbs = old_abi.keys()
    new_pdbs = new_abi.keys()
    overlap_pdbs = list(set(old_pdbs) & set(new_pdbs))
    print "old, new, overlap"
    print "{}, {}, {}".format(len(old_pdbs), len(new_pdbs), len(overlap_pdbs))

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

    return

if __name__ == "__main__":
    compare("/Users/jjeliazkov/Rosetta/tools/antibody-update/info/", "/Users/jjeliazkov/Rosetta/tools/antibody/info/")
