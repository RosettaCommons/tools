#!/usr/bin/env python

######################################################################
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options

if(is_release_mode()==False): #Not yet RELEASED!
	from SWA_dagman_python.utility.DAGMAN_util import setup_chemical_shift_args
######################################################################
from sys import argv,exit
import sys
import traceback
from glob import glob
import string
from os.path import basename, dirname, exists, expanduser
from time import sleep
from sets import Set
######################################################################


#extend_silent_struct_with_idealized_helix.py -extension_helix_sequence  ccgg -in_silent_file rebuild_bulge_region_FINAL.out

START_argv=copy.deepcopy(argv)

extension_helix_sequence=parse_options( argv, "extension_helix_sequence", "" )

if(extension_helix_sequence==""): error_exit_with_message("extension_helix_sequence==\"\"")

extension_helix_length=len(extension_helix_sequence)

if((extension_helix_length%2)!=0): error_exit_with_message("uneven extension_helix_length (%s) !" %(extension_helix_length))

strand_1_seq_num=range(1, (extension_helix_length/2)+1)
strand_2_seq_num=range((extension_helix_length/2)+1, extension_helix_length+1)

for base_ID in extension_helix_sequence:

	if(base_ID not in ["a", "g", "u", "c"]): error_exit_with_message("Error invalid base_ID (%s) in extension_helix_sequence (%s)" %(base_ID, extension_helix_sequence) )

in_silent_file=parse_options(argv, "in_silent_file", "" )

#out_silent_file=parse_optiona(argv, "out_silent_file", "" )

if(in_silent_file==""): error_exit_with_message("in_silent_file==\"\"")

if(exists(in_silent_file)==False): error_exit_with_message("in_silent_file (%s) doesn't exists!" %(in_silent_file))


if(in_silent_file[-4:]!=".out"): error_exit_with_message("in_silent_file[-4:]!=\".out\"")

out_silent_file="%s_with_extend_%s_helix.out" %(basename(in_silent_file[:-4]), extension_helix_sequence)

if(exists(out_silent_file)==True): 
	print "WARNING out_silent_file (%s) already exist .... removing!" %(out_silent_file)
	submit_subprocess("rm %s" %(out_silent_file))

##1. OK first create the idealized helix!##


extension_helix_pdb="extension_helix_%s.pdb" %(extension_helix_sequence)

create_extension_helix  ="create_idealize_helix_general.py -sequence %s -verbose True -VDW_screener_type NONE " %(extension_helix_sequence)
create_extension_helix +="-strand_1_seq_num %s -strand_2_seq_num %s " %( list_to_string(strand_1_seq_num), list_to_string(strand_2_seq_num) ) 
create_extension_helix +="-helix_filename %s " %(extension_helix_pdb)
submit_subprocess(create_extension_helix) 

if(exists(extension_helix_pdb)==False): error_exit_with_message("extension_helix_pdb (%s) doesn't exist!" %(extension_helix_pdb))

##2. Call Rosetta####

if(use_new_src_code()): 
	EXE = get_rosetta_EXE_specified_no_graphic_string("swa_rna_util", no_graphic_string)  #Feb 14, 2012: Haven't implement this for the updated code yet!
else:
	EXE = get_rosetta_EXE_specified_no_graphic_string("parin_test", no_graphic_string) 

database_folder= get_rosetta_database_folder() 

################################################################################################################



command=EXE
command += " -algorithm extend_silent_struct_with_idealized_helix"
command += " -database %s " %(database_folder)
command += " -in:file:silent %s " %(in_silent_file)
command += " -out:file:silent %s " %(out_silent_file)
command += " -in::file::tags %s " %(extension_helix_pdb)

print command 
submit_subprocess( command )

























print "-------------------------------------------------------------------------------------------"
print "Successfully RAN: %s" %(list_to_string(START_argv))
print "-------------------------------------------------------------------------------------------"

