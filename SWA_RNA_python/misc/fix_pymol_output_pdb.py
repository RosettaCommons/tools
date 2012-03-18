#!/usr/bin/python

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
import os
from os.path import basename, dirname, exists, expanduser
from time import sleep
import copy
import string
from os import popen 


from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options


#fix_pymol_output_pdb.py -input_pdb_file   3nj6_partner_RNA_A.pdb

#fix_pymol_output_pdb.py -input_pdb_file C4G5_aligned_lower_VDW_rep_screener.pdb

#AU_BP_aligned_to_AB_corner.pdb  
#A4U5_aligned_upper_VDW_rep_screener.pdb 
#C4G5_aligned_lower_VDW_rep_screener.pdb

input_pdb_file= parse_options( argv, "input_pdb_file", "" )

if(input_pdb_file==""): error_exit_with_message("input_pdb_file==\"\"")  

if(exists(input_pdb_file)==False): error_exit_with_message("input_pdb_file (%s) doesn't exist!!" %(input_pdb_file))  

output_pdb_file="FIXED_%s" %(input_pdb_file)


	
	

job_string_list=open(input_pdb_file).readlines()


if(exists(output_pdb_file)):
	print "Warning %s already exist...removing!" %(output_pdb_file)
	submit_subprocess( "rm %s " %(output_pdb_file) )

OUPUT_PDB = open( output_pdb_file, 'w')

for line in job_string_list:


	if(line=="END\n"): 
		#OUPUT_PDB.write("END\n")
		continue

	if(len(line)!= 81): error_exit_with_message("len(line)!= 81, len(line)=%d, line=%s" %(len(line), line) )  


	if(line[17]!='r'): error_exit_with_message("line[17]!='r', line=%s" %(line) )  
	if(line[18]!='G' and line[18]!='U' and line[18]!='C' and line[18]!='A'): error_exit_with_message("line[18]!='G' and line[18]!='U' and line[18]!='C' and line[18]!='A', line=%s" %(line))
	if(line[19]!=' '): error_exit_with_message("line[19]!=' ', line=%s" %(line) )  
	if(line[20]!=' '): error_exit_with_message("line[20]!=' ', line=%s" %(line) )  
	if(line[77]==' '): error_exit_with_message("line[77]==' ', line=%s" %(line) )

	new_line_list=[]

	for char in line:	
		new_line_list.append(char)
	
	new_line_list[17]=' '
	new_line_list[18]=line[17]
	new_line_list[19]=line[18]
	new_line_list[77]=' '

	seperator=""

	new_line=seperator.join(new_line_list)
	#print "new_line= %s" %(new_line)

	OUPUT_PDB.write(new_line)

OUPUT_PDB.close()

#consistency_check

rosetta_ready_pdb=output_pdb_file.lower().replace(".pdb", "_RNA_A.pdb")

if(output_pdb_file==rosetta_ready_pdb) : error_exit_with_message("output_pdb_file(%s)==rosetta_ready_pdb(%s)" %(output_pdb_file,rosetta_ready_pdb) )

if(exists(rosetta_ready_pdb)==True): error_exit_with_message("rosetta_ready_pdb (%s) already exist!!" %(rosetta_ready_pdb))  

submit_subprocess("SWA_make_rna_rosetta_ready.py %s" %(output_pdb_file))

if(exists(rosetta_ready_pdb)==False): error_exit_with_message("rosetta_ready_pdb (%s) doesn't exist!!" %(rosetta_ready_pdb))  

diff_line_list = popen('diff %s %s ' %(output_pdb_file, rosetta_ready_pdb) ).readlines()

if(len(diff_line_list)!=0):
	for diff_line in diff_line_list:
		print "diff_line: %s" %(diff_line)
	error_exit_with_message("len(diff_line_list)!=0")


submit_subprocess("rm %s " %(rosetta_ready_pdb) )


	#print line[17:20]

#convert line 1 to line 2
#line 1: 	ATOM     19  N6  rA  A   1      -1.054 -16.868  11.363  1.00  0.00           N 
#line 2: 	ATOM     19  N6   rA A   1      -1.054 -16.868  11.363  1.00  0.00
#		   		01234567890123456789012345678901234567890123456789012345678901234567890123456789	
#					0         1         2         3         4         5         6         7     
