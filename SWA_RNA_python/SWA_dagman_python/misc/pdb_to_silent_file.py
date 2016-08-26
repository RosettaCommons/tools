#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
import os
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################

#pdb_to_silent_file.py -pdb  s_1_RNA_A.pdb  s_2_RNA_A.pdb s_3_RNA_A.pdb s_4_RNA_A.pdb s_5_RNA_A.pdb s_6_RNA_A.pdb s_7_RNA_A.pdb s_8_RNA_A.pdb  s_9_RNA_A.pdb s_10_RNA_A.pdb  s_11_RNA_A.pdb s_12_RNA_A.pdb -output_silent_file  s_ALL.out 

#pdb_to_silent_file.py -pdb_glob_folder 4ee02ce090ce8/  -output_silent_file FR3D_All_tetraloop_3_angstrom.out

#pdb_to_silent_file.py -pdb rosetta_1gid_RNA_A_chain_A.pdb  -output_silent_file Virtualize_A24_1GID_A.out -virtual_res 22 -tag_name Virtual_A22_1GID_A

#pdb_to_silent_file.py -pdb S_000000.pdb  -output_silent_file Virtualize_U19_S_000000.out -virtual_res 12 19 -tag_name Virtual_U19_S_0

#pdb_to_silent_file.py -pdb S_000000.pdb  -output_silent_file Virtualize_U19_S_000000.out -virtual_res 12 19 -tag_name Virtual_U19_S_0

#pdb_to_silent_file.py -pdb 1_4_12_16_S_000000.pdb  -output_silent_file 1_4_12_16_S_000000.out -virtual_res 7
#pdb_to_silent_file.py -pdb 5_11_S_000000.pdb  -output_silent_file 5_11_S_000000.out -virtual_res 2 -virtual_ribose 1

#pdb_to_silent_file.py -pdb 3p59_RNA_A.pdb RENUM_DASLAB_TS002_1.pdb RENUM_DASLAB_TS002_2.pdb RENUM_DASLAB_TS002_3.pdb RENUM_DASLAB_TS002_4.pdb RENUM_DASLAB_TS002_5.pdb  -output_silent_file Nano_corner_SWA_models.out 

#pdb_to_silent_file.py -pdb  inverted_upper_daslab_ts001_1.pdb -tag_name inverted_upper_daslab_ts001_1 -output_silent_file inverted_upper_daslab_ts001_1.out -virtual_res 13

#pdb_to_silent_file.py -pdb  LOWER_daslab_ts001_3.pdb -tag_name LOWER_daslab_ts001_3 -output_silent_file LOWER_daslab_ts001_3.out -virtual_res 13

#pdb_to_silent_file.py -pdb  REGION_0_1_S_000000.pdb -tag_name S_0 -output_silent_file region_0_1_sample.cluster.out

#pdb_to_silent_file.py -pdb REGION_11_0_S_000000.pdb -tag_name S_0  -output_silent_file region_11_0_sample.cluster.out
######################################################################



no_graphic_string= parse_options( argv, "no_graphic", "" )



if(use_new_src_code()): 
	EXE = get_rosetta_EXE_specified_no_graphic_string("swa_rna_util", no_graphic_string) 
else:
	EXE = get_rosetta_EXE_specified_no_graphic_string("parin_test", no_graphic_string) 

database_folder= get_rosetta_database_folder() 

####################################################################

pdb_file_list =	  parse_options( argv, "pdb", [""] )
pdb_glob_folder = parse_options( argv, "pdb_glob_folder", "" )
list_of_virtual_res =  parse_options( argv, "list_of_virtual_res", [""] ) 
list_of_energy=  parse_options( argv, "list_of_energy", [""] ) 

tag_name =	  parse_options( argv, "tag_name", "" )
output_silent_file = parse_options( argv, "output_silent_file", "" )
virtual_res = parse_options( argv, "virtual_res", [-1] )
virtual_ribose = parse_options( argv, "virtual_ribose", [-1] )

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )

if(pdb_file_list!=[""]): 

	if(pdb_glob_folder!=""): error_exit_with_message("pdb_file_list!=[\"\"] and pdb_glob_folder!=\"\"")


elif(pdb_glob_folder!=""):

	if(exists(pdb_glob_folder)==False): error_exit_with_message("pdb_glob_folder (%s) doesn't exist!")

	pdb_file_list= glob("%s/*.pdb" %(pdb_glob_folder))
	pdb_file_list.sort()

else:
	error_exit_with_message('pdb_file_list==[""] and pdb_glob_folder==""')

for pdb_file in pdb_file_list:
	if(exists(pdb_file)==False): error_exit_with_message("pdb_file (%s) doesn't exist!" %(pdb_file) )

if(list_of_virtual_res!=[""]):
	if(len(pdb_file_list)!=len(list_of_virtual_res)): error_exit_with_message("len(pdb_file_list)!=len(list_of_virtual_res)")
	if(len(virtual_res)!=0): error_exit_with_message("len(list_of_virtual_res)>0 but len(virtual_res)!=0")

if(list_of_energy!=[""]):
	if(len(pdb_file_list)!=len(list_of_energy)): error_exit_with_message("len(pdb_file_list)!=len(list_of_energy)")



if(output_silent_file==""):
	if(len(pdb_file_list)==1):
		output_silent_file=basename(pdb_file_list[0]).replace(".pdb", ".out")
	else:
		error_exit_with_message("User need to pass in output_silent_file if len(pdb_file_list)!=1")

output_silent_file=os.path.abspath(output_silent_file)

if(exists(output_silent_file)):
	print "Warning output_silent_file already exist .. removing .." 
	submit_subprocess("rm %s " %(output_silent_file) )


#if( len(pdb_file.split()) != 1 ): error_exit_with_message("len(pdb_file.split()!=1")


####################################################################

command=EXE
command += " -database %s " %(get_rosetta_database_folder())
command += " -algorithm pdb_to_silent_file"
command += ' -s %s' %(list_to_string(pdb_file_list))
command += ' -output_silent_file %s ' %(output_silent_file)

if(len(virtual_res)!=0): command += ' -virtual_res %s ' %(list_to_string(virtual_res) )

if(list_of_virtual_res!=[""]): command += ' -list_of_virtual_res %s ' %(list_to_string(list_of_virtual_res) )

if(list_of_energy!=[""]): command += ' -list_of_energy %s ' %(list_to_string(list_of_energy) )

if(len(virtual_ribose)!=0): command += ' -virtual_ribose %s ' %(list_to_string(virtual_ribose) )

if(tag_name!=""):
	if(len(pdb_file_list)!=1): error_exit_with_message('tag_name!="" and len(pdb_file_list)!=1' )
	command += ' -tag_name %s' %(tag_name)

command += ' > pdb_to_silent_file_output.txt '

print command 
submit_subprocess( command )


if(exists(output_silent_file)==False): error_exit_with_message("output_silent_file (%s) doesn't exist!" %(output_silent_file) )

#if(output_silent_file!=""):
#	default_silent_file=pdb_file.replace(".pdb", ".out")
#	if(exists(default_silent_file)==False): error_exit_with_message("default_silent_file (%s) doesn't exist!" %(default_silent_file) )
#	submit_subprocess( "mv %s %s" %(default_silent_file, output_silent_file) )
