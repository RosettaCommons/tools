#!/usr/bin/env python

from os import system, chdir, getcwd, path
from sys import argv, exit, stderr
from os.path import exists, isdir


def get_structure_profile_checkpoint( template_fn ):

    ## get PDB2VALL_PATH
    PDB2VALL_PATH = path.abspath(path.dirname(__file__))
    PDB2VALL_PATH = PDB2VALL_PATH[0:PDB2VALL_PATH.rfind("/pdb2vall/")];
    PDB2VALL_PATH += "/pdb2vall/"
    if not PDB2VALL_PATH:
        stderr.write("ERROR: you should specify the path where your packages are first.\n")
        return 0

    # ALL THE SCRIPTS_FN NECESSARY
    script_get_depth_data                         = PDB2VALL_PATH + "structure_profile_scripts/make_depthfile.py"
    script_make_sequence_fragments                = PDB2VALL_PATH + "structure_profile_scripts/make_sequence_fragments.pl"
    script_get_single_chain_pdb                   = PDB2VALL_PATH + "pdb_scripts/get_pdb_new.py"
    script_get_fasta_from_pdb                     = PDB2VALL_PATH + "pdb_scripts/pdb2fasta.py"
    script_make_alignment_from_fragfile           = PDB2VALL_PATH + "structure_profile_scripts/make_alignment_from_fragfile.pl"
    script_create_checkpoint_from_fasta_alignment = PDB2VALL_PATH + "structure_profile_scripts/create_checkpoint_from_fasta_alignment.pl"

    base_dir = getcwd()
    working_dir = getcwd() + "/structure_profile_checkpoint/"

    # MKCD A DIRECTORY - STRUCTURE_PROFILE_CHECKPOINT/
    if not isdir("structure_profile_checkpoint"):
        system("mkdir structure_profile_checkpoint/")

    chdir( working_dir )

    depth_results = template_fn[:5] + ".depth-residue.depth"
    if not exists( depth_results ):
        cmd = script_get_depth_data + " " + template_fn
        print cmd + "\n"
        system( cmd )

    cmd = script_get_single_chain_pdb + " " + template_fn[:4] + " " + template_fn[4]
    print cmd + "\n"
    system( cmd )

    fasta_fn = "%s.fasta" % template_fn[:5]
    cmd = script_get_fasta_from_pdb + " " + template_fn + " > %s.fasta" % template_fn
    print cmd + "\n"
    system( cmd )

    # THE DISTINGUISH BETWEEN THE FASTA FROM PDB_SEQRES.TXT AND THE ONE FROM PDB2FASTA.PY
    if (exists(fasta_fn)):
        system("rm %s" % ( fasta_fn ))
    system("ln -s %s.fasta %s" % ( template_fn, fasta_fn ))

    sequence_fragments     = template_fn[:5] + ".50.9mers"
    sequence_fragments_fsc = template_fn[:5] + "_frags.fsc.50.9mers"

    if not exists( sequence_fragments and sequence_fragments_fsc ):
        print "get_structure_profile_checkpoint(): running fragment_picker.app to pick sequence fragments"

        cmd = script_make_sequence_fragments + " -verbose -n_frags 50 -n_candidates 500 -frag_sizes 9 -depth %s -native %s.pdb %s.fasta" % ( depth_results, template_fn[:5], template_fn[:5] )
        print cmd + "\n"
        system( cmd )

    if not exists( sequence_fragments and sequence_fragments_fsc ):
        stderr.write("get_structure_profile_checkpoint(): something wrong with the the make_sequence_fragments script\n"); exit()

    if not exists( sequence_fragments + ".ali.fasta" ):
        cmd = script_make_alignment_from_fragfile + " " + fasta_fn
        print cmd + "\n"
        system( cmd )

    if not exists( sequence_fragments + ".ali.fasta.new.blast.checkpoint" ):
        cmd = script_create_checkpoint_from_fasta_alignment + " " + fasta_fn + " " + sequence_fragments + ".ali.fasta"
        print cmd + "\n"
        system( cmd )

    chdir( base_dir )

    cmd = "ln -s ./structure_profile_checkpoint/" + sequence_fragments + ".ali.fasta.new.blast.checkpoint ."
    print cmd + "\n"
    system( cmd )

if __name__ == "__main__":
    if len( argv ) < 2:
        print
        print "USAGE: %s <template_id: eg. 2oxgZ.pdb>" % argv[0]
        print
        print "-"*75
        exit()

    print get_structure_profile_checkpoint( argv[1] )

