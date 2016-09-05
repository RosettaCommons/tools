#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
from glob import glob
import os
import copy

from SWA_dagman_python.utility.PDB_operations import *
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options

##############################USAGE#################################################################

#create_torsion_database_check_length.py -PDB_src_table PDB_src_table.txt

####################################################################################################


copy_argv=copy.deepcopy(argv)

#####################################
print "Enter %s " %(list_to_string(copy_argv) )

HOMEDIR = expanduser('~') 

print "HOMEDIR= %s" %(HOMEDIR)

#############################################

PDB_src_table=parse_options( argv, "PDB_src_table", "")

if(exists(PDB_src_table)==False): error_exit_with_message("PDB_src_table (%s) doesn't exist!" %(PDB_src_table))

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )

#####################################################################

PDB_SRC_LINES=PDB_LINES = open( PDB_src_table ).readlines()

first_line=True

total_nts=0
num_RNA_struct=0

for pdb_src_line in PDB_SRC_LINES:

	if(first_line):
		first_line=False
		continue

	pdb_file=pdb_src_line.split()[0]

	if(exists(pdb_file)==False): error_exit_with_message("pdb_file (%s) doesn't exist!" %(pdb_file))

	pdb_length=int(pdb_src_line.split()[1])

	full_length_torsion_outfile="full_length_rosetta_rs_ready_" + pdb_file.replace(".pdb", "") + ".torsions"

	if(exists(full_length_torsion_outfile)==False): error_exit_with_message("full_length_torsion_outfile (%s) doesn't exist!" %(full_length_torsion_outfile))

	torsion_file_lines=open( full_length_torsion_outfile ).readlines()

	if(pdb_length!=len(torsion_file_lines)): error_exit_with_message("pdb_length!=len(torsion_file_lines) for pdb_file (%s)" %(pdb_file))

	total_nts+=pdb_length
	num_RNA_struct+=1

summary_line= "SUMMARY: NUM_RNA_STRUCT=%s | TOTAL_NTS=%s" %(num_RNA_struct, total_nts)

print summary_line

