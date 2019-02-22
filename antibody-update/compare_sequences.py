"""
Compare sequences from grafted models to input FASTA.
"""
import os
import sys
import pyrosetta

if __name__ == "__main__":
    pyrosetta.init("-mute all")

    # target directory (containing dirs containing models)
    # e.g. target_dir/1abc/grafting, target_dir/1def/grafting, ...
    # and output file for analysis
    target_dir = "benchmark_abs_old_db"
    #outfile = "old_ab_seq_diffs.csv"

    # read targets by PDB ID
    os.chdir("/Users/lqtza/Rosetta/tools/antibody-update")
    with open("pdbids.txt", "r") as inf:
        pdbs = [x.strip() for x in inf.readlines()]

    # for each PDB ID, read the FASTA and model-0
    for pdb in pdbs:

        print(pdb)

        # skip some pdbs not in PDB
        # not in new db: 3umt, 3nps, 4h0h
        #if pdb in ["1x9q", "3eo9", "3ifl"]: continue
        if pdb in ["1mfa"]: continue

        # load crystal from database (equivalent to fASTA, I hope)
        crystal = "benchmark_crystals/{}_trunc.pdb".format(pdb)
        crystal_pose = pyrosetta.pose_from_file(crystal)
        crystal_seq = crystal_pose.sequence()

        # load models and compare
        model = "{}/{}/grafting/model-0.pdb".format(target_dir, pdb)
        model_pose = pyrosetta.pose_from_file(model)
        model_seq = model_pose.sequence()

        # compare sequences
        if crystal_seq != model_seq:
            print "xtal:  {}".format(crystal_seq)
            print "model: {}".format(model_seq)

        #for res in range(crystal_pose.size()):
            # convert to pdb numbering and check if both exist
