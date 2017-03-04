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


#Rosetta_import_and_dump_pdb.py  -s 1gid_RNA_A.pdb -no_graphic True

print "Enter: " , list_to_string(argv)

pose_name= parse_options( argv, "s", "" )

if(pose_name==""): error_exit_with_message('pose_name==""')

no_graphic_string= parse_options( argv, "no_graphic", "" )

reset_o2star_torsion= parse_options( argv, "reset_o2star_torsion", "False" )

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )

####################################################################



if(use_new_src_code()): 
	EXE = get_rosetta_EXE_specified_no_graphic_string("swa_rna_util", no_graphic_string) 
else:
	EXE = get_rosetta_EXE_specified_no_graphic_string("parin_test", no_graphic_string) 


database_folder= get_rosetta_database_folder() 

#####################################################################


print "input_pdb=%s" %(pose_name)

command=EXE
command += " -algorithm import_and_dump_pdb"
command += " -database %s " %(database_folder)
command += " -s " + pose_name

print command 
submit_subprocess( command )




