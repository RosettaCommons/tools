#!/usr/bin/env python

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


#remove_hydrogen_from_pdb.py -s  S_000004.pdb

input_pdb_file= parse_options( argv, "s", "" )

output_pdb_file= parse_options( argv, "output_pdb", "")



if(input_pdb_file==""): error_exit_with_message("input_pdb_file==\"\"")  

if(exists(input_pdb_file)==False): error_exit_with_message("input_pdb_file (%s) doesn't exist!!" %(input_pdb_file))  


######################################################################################################
if(output_pdb_file==""): 
	output_pdb_file="remove_hydrogen_%s" %(basename(input_pdb_file))
else:
	"User inputted output_pdb_file= %s " %(output_pdb_file)
	if(exists(output_pdb_file)): error_exit_with_message("User input output_pdb_file (%s) already exist!" %(output_pdb_file) )

if(exists(output_pdb_file)):
	print "Warning %s already exist...removing!" %(output_pdb_file)
	submit_subprocess( "rm %s " %(output_pdb_file) )

OUPUT_PDB = open( output_pdb_file, 'w')
######################################################################################################

atom_count=0

job_string_list=open(input_pdb_file).readlines()

for line in job_string_list:

	

	if(line=="END\n"): continue 

	if(line[0:3]=="TER"): continue 

		
	#if(len(line)!= 81): error_exit_with_message("len(line)!= 81, len(line)=%d, line=%s" %(len(line), line) )  

	if(line[18]!='r'): error_exit_with_message("line[18]!='r', line=%s" %(line) )  
	if(line[19]!='G' and line[19]!='U' and line[19]!='C' and line[19]!='A'): error_exit_with_message("line[19]!='G' and line[19]!='U' and line[19]!='C' and line[19]!='A' line=%s" %(line))
	if(line[20]!=' '): error_exit_with_message("line[20]!=' ', line=%s" %(line) )  

	new_line_list=[]

	for char in line:	new_line_list.append(char)


	if(line[12:16]=="1H5*"): continue
	if(line[12:16]=="2H5*"): continue
	if(line[12:16]==" H4*"): continue
	if(line[12:16]==" H3*"): continue
	if(line[12:16]==" H1*"): continue
	if(line[12:16]=="1H2*"): continue
	if(line[12:16]=="2HO*"): continue

	if(line[12:16]=="1H2 "): continue
	if(line[12:16]=="2H2 "): continue
	if(line[12:16]=="1H4 "): continue
	if(line[12:16]=="2H4 "): continue
	if(line[12:16]=="1H6 "): continue
	if(line[12:16]=="2H6 "): continue

	if(line[12:16]==" H1 "): continue
	if(line[12:16]==" H3 "): continue
	if(line[12:16]==" H2 "): continue
	if(line[12:16]==" H5 "): continue
	if(line[12:16]==" H6 "): continue
	if(line[12:16]==" H8 "): continue


	if(line.count("H")!=0): error_exit_with_message("line.count(\"H\")!=0, line=%s" %(line))

	atom_count+=1

	new_line_list[6:11]="%5s" %(atom_count)


	seperator=""

	new_line=seperator.join(new_line_list)
	#print "new_line= %s" %(new_line)

	OUPUT_PDB.write(new_line)

OUPUT_PDB.close()

'''
ATOM   1575  O2P  rA A  74      36.620  69.335  42.589  1.00 10.00           O  
ATOM     24 1H5*  rG R   1      98.817 -44.987  15.431  1.00  0.00              
ATOM     25 2H5*  rG R   1     100.046 -44.552  14.237  1.00  0.00              
ATOM     26  H4*  rG R   1      97.690 -44.283  13.549  1.00  0.00              
ATOM     27  H3*  rG R   1      99.544 -41.837  13.360  1.00  0.00              
ATOM     28  H1*  rG R   1      95.744 -41.613  13.960  1.00  0.00              
ATOM     29 1H2*  rG R   1      97.656 -40.736  12.212  1.00  0.00              
ATOM     30 2HO*  rG R   1      96.198 -41.693  11.150  1.00  0.00              
ATOM     31  H1   rG R   1      95.469 -35.679  13.910  1.00  0.00              
ATOM     32 1H2   rG R   1      93.818 -37.967  11.925  1.00  0.00              
ATOM     33 2H2   rG R   1      94.022 -36.309  12.381  1.00  0.00              
ATOM     34  H8   rG R   1      98.782 -40.770  16.074  1.00  0.00 
01234567890123456789012345678901234567890123456789012345678901234567890123456789	
'''
#					ATOM      2  C5*  rG A   1      40.311  50.542  48.455  1.00 33.79           C  (Rosetta)

#		   		01234567890123456789012345678901234567890123456789012345678901234567890123456789	
#					0         1         2         3         4         5         6         7     
