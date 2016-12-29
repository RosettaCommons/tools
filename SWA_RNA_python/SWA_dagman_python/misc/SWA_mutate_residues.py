#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################

output_pdb= parse_options( argv, "output_pdb", "" )
no_graphic_string= parse_options( argv, "no_graphic", "" )
verbose= parse_options( argv, "verbose", "True" )

if(len(argv)!=3): error_exit_with_message("len(argv)!=3, leftover_argv=%s" %(list_to_string(argv) ) )

pose_name = argv[1]
mutate_residue_string= argv[2]

if(verbose): print pose_name

if(exists(pose_name)==False): error_exit_with_message("pose_name (%s) doesn't exist!" %(pose_name) )

##############################################################

if(use_new_src_code()):
	command=get_rosetta_EXE_specified_no_graphic_string("swa_rna_util", no_graphic_string)  
	command += " -algorithm mutate_residues"
else:
	command=get_rosetta_EXE_specified_no_graphic_string("parin_test", no_graphic_string)  
	command += " -algorithm Mutate_residue"

command += " -rebuild_sequence " + mutate_residue_string
command += " -database %s" %(get_rosetta_database_folder())
command += " -s " + pose_name
command += " > mutate_" + pose_name + "_" + mutate_residue_string + "_output.txt"
if(verbose): print command 
popen( command )

#mutate_G1-G2_expand_radius_100_1s72_RNA_A_35_40_3_8.pdb
rosetta_default_output_pdb="mutate_%s_%s" %(mutate_residue_string, pose_name)

if(exists(rosetta_default_output_pdb)==False): error_exit_with_message("rosetta_default_output_pdb (%s) doesn't exist!" %(rosetta_default_output_pdb) ) 

if(output_pdb!=""): 
	print "mv %s %s " %(rosetta_default_output_pdb, output_pdb)
	submit_subprocess("mv %s %s " %(rosetta_default_output_pdb, output_pdb) )
