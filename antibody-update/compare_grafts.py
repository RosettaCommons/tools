"""
Compare grafted models to crystal structures.
To be used in conjunction with the grafting benchmarks.
A single run will compare grafted models to crystal structures.
Output location should be specified
"""
import os
import sys

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

def get_framework_residue_maps(crystal_pose, model_pose):
    # atom maps for alignment
    # we want to map from crystal to model, then align crystal to model
    # define conserved framework regions
    #conserved_frh_residues = [10,11,12,13,14,15,16,17,18,19,20,21,21,23,24,25,36,37,38,39,40,41,42,43,44,45,46,47,48,49,66,69,70,71,72,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,103,104,105]
    conserved_frh_residues = [10,11,12,13,14,15,16,17,18,19,20,21,21,23,24,25,36,37,38,39,40,42,43,44,45,46,47,48,49,66,69,70,71,72,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,103,104,105]
    #conserved_frl_residues = [10,11,12,13,14,15,16,17,18,19,20,21,21,23,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,57,58,59,60,61,62,63,64,65,66,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,98,99,100]
    conserved_frl_residues = [11,12,13,14,15,16,17,18,19,20,21,21,23,35,36,37,38,39,41,42,43,44,45,46,47,48,49,57,58,59,60,61,62,63,64,65,66,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,98,99,100]

    frl_map = pyrosetta.rosetta.core.id.AtomID_Map_AtomID()
    pyrosetta.rosetta.core.pose.initialize_atomid_map(frl_map, crystal_pose, pyrosetta.rosetta.core.id.AtomID.BOGUS_ATOM_ID())

    frh_map = pyrosetta.rosetta.core.id.AtomID_Map_AtomID()
    pyrosetta.rosetta.core.pose.initialize_atomid_map(frh_map, crystal_pose, pyrosetta.rosetta.core.id.AtomID.BOGUS_ATOM_ID())

    # convert conserved numbering to pose numbering
    # there should be no 0s but check anyway
    for res in conserved_frh_residues:

        cres = crystal_pose.pdb_info().pdb2pose("H", res)
        mres = model_pose.pdb_info().pdb2pose("H", res)

        if cres == 0 or mres == 0:
            sys.exit("Missing conserved FRH residue {}!".format(res))

        # otherwise map backbone of crystal to model
        # first four atoms are N, CA, C, O
        for i in [1,2,3,4]:
            cid = pyrosetta.rosetta.core.id.AtomID(i, cres)
            mid = pyrosetta.rosetta.core.id.AtomID(i, mres)
            frh_map.set(cid, mid)

    # convert conserved numbering to pose numbering
    # there should be no 0s but check anyway
    for res in conserved_frl_residues:

        cres = crystal_pose.pdb_info().pdb2pose("L", res)
        mres = model_pose.pdb_info().pdb2pose("L", res)

        if cres == 0 or mres == 0:
            sys.exit("Missing conserved FRL residue {}!".format(res))

        # otherwise map backbone of crystal to model
        # first four atoms are N, CA, C, O
        for i in [1,2,3,4]:
            cid = pyrosetta.rosetta.core.id.AtomID(i, cres)
            mid = pyrosetta.rosetta.core.id.AtomID(i, mres)
            frl_map.set(cid, mid)

    return frl_map, frh_map

def get_cdr_atoms(crystal_pose, model_pose, crystal_abi, model_abi):
    # all alignment in rosetta is a headache -- lets just pull out the residues
    # we care about and calc rms on our own

    crystal_cdr_dict = {}
    model_cdr_dict = {}

    for enum in [antibody.h1, antibody.h2, antibody.h3, antibody.l1, antibody.l2, antibody.l3]:
        # init vecotrs to store atom coords
        crystal_cdr_dict[enum] = pyrosetta.rosetta.utility.vector1_numeric_xyzVector_double_t()
        model_cdr_dict[enum] = pyrosetta.rosetta.utility.vector1_numeric_xyzVector_double_t()

        # get CDR residues and check for length
        cloop = crystal_abi.get_CDR_loop(enum, crystal_pose)
        mloop = model_abi.get_CDR_loop(enum, model_pose)

        if cloop.size() != mloop.size():
            sys.exit("Cannot compare CDRs of different length")

        # store paired coords
        for i in range(cloop.size()+1): # we should include the end

            # first four atoms are N, CA, C, O
            for j in [1,2,3,4]:

                cid = pyrosetta.rosetta.core.id.AtomID(j, cloop.start()+i)
                mid = pyrosetta.rosetta.core.id.AtomID(j, mloop.start()+i)

                crystal_cdr_dict[enum].append(crystal_pose.xyz(cid))
                model_cdr_dict[enum].append(model_pose.xyz(mid))

    return crystal_cdr_dict, model_cdr_dict

def compare_pair(crystal_pose, crystal_abi, model_pose, model_abi):
    rmsds = {"frh":0.0, "frl":0.0, "l1":0.0, "l2":0.0, "l3":0.0, "h1":0.0, "h2":0.0, "h3":0.0, "ocd":0.0}

    # get OCD
    crystalO = antibody.vl_vh_orientation_coords(crystal_pose, crystal_abi)
    modelO = antibody.vl_vh_orientation_coords(model_pose, model_abi)
    rmsds["ocd"] = calc_orientation_diff(crystalO, modelO)

    # get maps for rmsd calc
    frl_map, frh_map = get_framework_residue_maps(crystal_pose, model_pose)
    crystal_cdr_dict, model_cdr_dict = get_cdr_atoms(crystal_pose, model_pose, crystal_abi, model_abi)

    # align and calc rmsd for each chain
    # function uses map from first arg to second arg (double check?)
    rmsds["frh"] = scoring.superimpose_pose(crystal_pose, model_pose, frh_map)

    # get cdr rmsds while we're at it
    rmsds["h1"] = pyrosetta.rosetta.numeric.model_quality.calc_rms(crystal_cdr_dict[antibody.h1], model_cdr_dict[antibody.h1])
    rmsds["h2"] = pyrosetta.rosetta.numeric.model_quality.calc_rms(crystal_cdr_dict[antibody.h2], model_cdr_dict[antibody.h2])
    rmsds["h3"] = pyrosetta.rosetta.numeric.model_quality.calc_rms(crystal_cdr_dict[antibody.h3], model_cdr_dict[antibody.h3])

    # repeat for L -- maybe make this a function too...
    rmsds["frl"] = scoring.superimpose_pose(crystal_pose, model_pose, frl_map)
    # get cdr rmsds while we're at it
    rmsds["l1"] = pyrosetta.rosetta.numeric.model_quality.calc_rms(crystal_cdr_dict[antibody.l1], model_cdr_dict[antibody.l1])
    rmsds["l2"] = pyrosetta.rosetta.numeric.model_quality.calc_rms(crystal_cdr_dict[antibody.l2], model_cdr_dict[antibody.l2])
    rmsds["l3"] = pyrosetta.rosetta.numeric.model_quality.calc_rms(crystal_cdr_dict[antibody.l3], model_cdr_dict[antibody.l3])

    return rmsds

if __name__ == "__main__":
    # imports and inits
    import pyrosetta
    from pyrosetta.rosetta.protocols import antibody
    from pyrosetta.rosetta.core import scoring
    pyrosetta.init()

    # target directory (containing dirs containing models)
    # e.g. target_dir/1abc/grafting, target_dir/1def/grafting, ...
    # and output file for analysis
    target_dir = "benchmark_abs_old_db"
    outfile = "test.csv"

    # number of models (typically 10)
    n_models = 10

    # read targets by PDB ID
    os.chdir("/Users/lqtza/Rosetta/tools/antibody-update")
    with open("pdbids.txt", "r") as inf:
        pdbs = [x.strip() for x in inf.readlines()]

    outf = open(outfile, "w")
    outf.write("pdb, model, ocd, frh, frl, h1, h2, h3, l1, l2, l3\n")
    for pdb in pdbs:
        print(pdb)
        # skip some pdbs not in PDB
        # not in new db: 3umt, 3nps, 4h0h
        if pdb in ["1mfa", "1x9q", "3eo9", "3ifl", "3mlr", "3nps"]: continue
        # load crystal from database
        crystal = "benchmark_crystals/{}_trunc.pdb".format(pdb)
        crystal_pose = pyrosetta.pose_from_file(crystal)
        crystal_abi = antibody.AntibodyInfo(crystal_pose)
        # load models and compare
        for i in range(n_models):
            model = "{}/{}/grafting/model-{}.pdb".format(target_dir, pdb, i)
            model_pose = pyrosetta.pose_from_file(model)
            model_abi = antibody.AntibodyInfo(model_pose)
            res = compare_pair(crystal_pose, crystal_abi, model_pose, model_abi)
            outf.write("{}, {}, ".format(pdb, i))
            outf.write("{}, ".format(res["ocd"]))
            outf.write("{}, ".format(res["frh"]))
            outf.write("{}, ".format(res["frl"]))
            outf.write("{}, ".format(res["h1"]))
            outf.write("{}, ".format(res["h2"]))
            outf.write("{}, ".format(res["h3"]))
            outf.write("{}, ".format(res["l1"]))
            outf.write("{}, ".format(res["l2"]))
            outf.write("{}\n".format(res["l3"]))
    outf.close()
