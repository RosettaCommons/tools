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
import time

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options


#remove_columns.py -remove_cols 3 11 22 23 24 25 -infile full_length_Rosetta_1jj2_NO_5S.torsions
#remove_columns.py -remove_cols 3 11 22 23 24 25 -infile 1jj2.torsions


remove_cols= parse_options( argv, "remove_cols", [0] ) 

if(len(remove_cols)==0): error_exit_with_message("len(remove_cols)==0")

infile=parse_options( argv, "infile", "" ) 

if(infile==""): error_exit_with_message("infile==\"\"")
if(exists(infile)==False): error_exit_with_message("infile (%s) doesn't exist!" %(infile) )

try:
	DATA = open(infile,'r')
except:
	error_exit_with_message("cannot open %s " %(infile) )

outfile="remove_cols%s_%s" %(list_to_string(remove_cols, "_"), basename(infile) ) 

OUTFILE=open(outfile, 'w')

while(True):
	line = DATA.readline()

	if(line==''): break #End of file!

	line_split=line.split()
	new_line=""

	for n in range(len(line_split)):

		curr_col=n+1

		if(curr_col in remove_cols): continue

		new_line+="   %10s" %(line_split[n])

	OUTFILE.write("%s\n" %(new_line))

DATA.close()
OUTFILE.close()
		
