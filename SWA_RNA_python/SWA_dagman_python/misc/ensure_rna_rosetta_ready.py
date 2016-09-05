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


#ensure_rna_rosetta_ready.py -input_pdb target1_RNA_A.pdb

input_pdb= parse_options( argv, "input_pdb", "" )

if(input_pdb==""): error_exit_with_message('input_pdb==""')

if(exists(input_pdb)==False): error_exit_with_message("input_pdb %s doesn't exist!" %(input_pdb) )

input_pdb=os.path.abspath(input_pdb)

temp_folder="TEMP_FOLDER/"

base_folder=os.path.abspath(".")

ID=1
while(exists(temp_folder)):
	temp_folder="TEMP_FOLDER_%d" %(ID)
	ID+=1

submit_subprocess("mkdir %s " %(temp_folder) )

os.chdir( "%s/" %(temp_folder) )

submit_subprocess("cp -r %s %s" %(input_pdb, basename(input_pdb) ) )

input_pdb=basename(input_pdb)

rosetta_ready_pdb="ROSETTA_READY_TEST_%s" %(input_pdb)

if(exists(rosetta_ready_pdb)==True): error_exit_with_message("rosetta_ready_pdb (%s) already exist!!" %(rosetta_ready_pdb))  

submit_subprocess("SWA_make_rna_rosetta_ready.py %s -output_pdb %s " %(input_pdb, rosetta_ready_pdb) )

if(exists(rosetta_ready_pdb)==False): error_exit_with_message("rosetta_ready_pdb (%s) doesn't exist!!" %(rosetta_ready_pdb))  

diff_line_list = popen('diff %s %s ' %(input_pdb, rosetta_ready_pdb) ).readlines()

if(len(diff_line_list)!=0):
	for diff_line in diff_line_list:
		print "diff_line: %s" %(diff_line)
	print "current_folder= %s " %(os.path.abspath("."))
	error_exit_with_message("len(diff_line_list)!=0")


os.chdir( "%s/" %(base_folder) )

submit_subprocess("rm -r %s " %(temp_folder) )
