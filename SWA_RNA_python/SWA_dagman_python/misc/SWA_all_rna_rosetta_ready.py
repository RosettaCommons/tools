#!/usr/bin/env python

import string
from sys import argv,stdout, stderr
from os import popen,system
from os.path import exists,dirname,basename,abspath
######################################################################

from SWA_dagman_python.database.SWA_amino_acids import longer_names
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################

##SWA_make_rna_rosetta_ready_wrapper.py -pdb_glob_folder 4ee02ce090ce8/ -output_folder RR_FORMATTED/

pdb_file_list =	  parse_options( argv, "pdb", [""] )
pdb_glob_folder = parse_options( argv, "pdb_glob_folder", "" )

output_folder=parse_options( argv, "output_folder", "" )

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )

if(pdb_file_list!=[""]): 

	if(pdb_glob_folder!=""): error_exit_with_message("pdb_file_list!=[\"\"] and pdb_glob_folder!=\"\"")


elif(pdb_glob_folder!=""):

	if(exists(pdb_glob_folder)==False): error_exit_with_message("pdb_glob_folder (%s) doesn't exist!")

	pdb_file_list= glob("%s/*.pdb" %(pdb_glob_folder))
	pdb_file_list.sort()

else:
	error_exit_with_message('pdb_file_list==[""] and pdb_glob_folder==""')

output_pdb= parse_options( argv, "output_pdb", "" )


alter_conform= parse_options( argv, "alter_conform", "A" )

if(alter_conform not in ['A', 'B']): error_exit_with_message("alter_conform not in ['A', 'B'] | alter_conform=%s" %(alter_conform))

if(output_folder!=""): 
	if(exists(output_folder)): 
		print "Warning output_folder (%s) already exist! ...removing..." %(output_folder)
		submit_subprocess("rm -r %s" %(output_folder))

	submit_subprocess("mkdir %s" %(output_folder))

for pdb in pdb_file_list:

	output_pdb="RR_%s" %(basename(pdb))

	command="SWA_make_rna_rosetta_ready.py %s -output_pdb %s -alter_conform %s " %(pdb, output_pdb, alter_conform)

	submit_subprocess(command)

	if(exists(output_pdb)==False): error_exit_with_message("output_pdb (%s) doesn't exist!" %(output_pdb))

	#submit_subprocess("Rosetta_import_and_dump_pdb.py  -s %s -no_graphic True" %(output_pdb))

	#Rosetta_import_and_dump_pdb="rosetta_%s" %(output_pdb)

	#if(exists(Rosetta_import_and_dump_pdb)==False): error_exit_with_message("Rosetta_import_and_dump_pdb (%s) doesn't exist!" %(Rosetta_import_and_dump_pdb))

	#submit_subprocess("mv %s %s " %(Rosetta_import_and_dump_pdb,  output_pdb))

	if(output_folder!=""): submit_subprocess("mv %s %s/%s " %(output_pdb, output_folder, output_pdb))





