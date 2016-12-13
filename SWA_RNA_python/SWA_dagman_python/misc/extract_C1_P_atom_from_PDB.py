#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options

#This code assume that input PDB is rna_rosetta_ready!

#extract_C1_P_atom_from_PDB.py -pdb 1q9a_RNA_A.pdb

pdb_file= parse_options( argv, "pdb", "")
include_first_N= parse_options( argv, "include_first_N", "True")

prefix=""
if(include_first_N):
	prefix="N_"

output_file= prefix + "C1_P_atom_"+pdb_file

pdb_string_list=open(pdb_file).readlines()
OUPUT_PDB = open( output_file, 'w')


for line in pdb_string_list:

	if(line=="END\n"): 
		#OUPUT_PDB.write("END\n")
		continue

	if(len(line)< 81): error_exit_with_message("len(line)< 81, len(line)=%d, line=%s" %(len(line), line) )  

	if(line[0:4] != 'ATOM'): error_exit_with_message("line[0:4] != 'ATOM', line=%s" %(line) ) 

	if(line[16]!=' '): error_exit_with_message("line[16]!=' ' != 'ATOM', line=%s" %(line) ) 

	if(line[13:16]=="C1*"): 
		line=line[0:13] + "C1'" + line[16:]

	write_line=False

	if(line[13:16]=="C1'" or line[13:16]=="P  "):
		write_line=True

	if(include_first_N):
		if(line[19]=='G' or line[19]=='A'):
			if(line[13:16]=="N9 "): write_line=True

		if(line[19]=='C' or line[19]=='U'):
			if(line[13:16]=="N1 "): write_line=True

	if(write_line):	
		#print line
		OUPUT_PDB.write(line)


#ATOM      2  C5*  rG A   1       3.726   3.822  22.615  1.00 84.57           C  
#0123456789012345678901234567890123456789012345678901234567890123456789
#0					1					2				 3					 4					5					6
