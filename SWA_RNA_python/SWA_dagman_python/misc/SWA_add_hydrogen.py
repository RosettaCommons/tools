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

#SWA_add_hydrogen.py USER_POSE_NAME -reset_o2star_torsion True

no_graphic_string= parse_options( argv, "no_graphic", "" )

reset_o2star_torsion= parse_options( argv, "reset_o2star_torsion", "False" )

if(len(argv)!=2): error_exit_with_message("len(argv)!=2, leftover_argv=%s" %(list_to_string(argv) ) )

pose_name = argv[1]

####################################################################



if(use_new_src_code()): 
	EXE = get_rosetta_EXE_specified_no_graphic_string("swa_rna_util", no_graphic_string) 
else:
	EXE = get_rosetta_EXE_specified_no_graphic_string("parin_test", no_graphic_string) 

database_folder= get_rosetta_database_folder() 

#####################################################################

print "input_pdb=%s" %(pose_name)


command=EXE
command += " -algorithm o2star_packer"
command += " -database %s " %(database_folder)
command += " -s " + pose_name
if(reset_o2star_torsion):
	command += " -reset_o2star_torsion true " 

command += " -constant_seed true -jran 1111111 " #Added this on Jan 24, 2012..This prevents randomness in the o2star_trail!

#command += " > add_hydrogen_" + pose_name + "_output.txt"
print command 
submit_subprocess( command )




