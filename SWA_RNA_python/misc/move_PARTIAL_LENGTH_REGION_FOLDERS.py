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

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
from SWA_dagman_python.parser.SWA_parse_internal_arguments import *

final_cluster_dag_file= parse_options( argv, "final_cluster_dag_file", "CONDOR/REGION_FINAL_cluster.condor")

full_length_FOLDER_list=get_full_length_REGION_FOLDER_list(final_cluster_dag_file, verbose=True)

folder_globstring = 'REGION_*_*'
all_FOLDER_list = glob( folder_globstring )
all_FOLDER_list.sort()


print
print "all_FOLDER_list=", all_FOLDER_list
print


partial_length_FOLDER_list=[]

#consistency check:
for FOLDER in all_FOLDER_list:
	FOLDER_split=FOLDER.split("_")
	if(len(FOLDER_split)!=3): error_exit_with_message("len(FOLDER_split)!=3, FOLDER_split=%s" %(list_to_string(FOLDER_split)) )

	
	#try:
	#	int(FOLDER_split[0])
	#except:
	#	error_exit_with_message("cannot convert FOLDER_split[0] to int, FOLDER_split[0]= %s, FOLDER= %s" %(FOLDER_split[0], FOLDER) )

	try:
		int(FOLDER_split[1])
	except:
		error_exit_with_message("cannot convert FOLDER_split[1] to int, FOLDER_split[1]= %s, FOLDER= %s" %(FOLDER_split[1], FOLDER) )

	try:
		int(FOLDER_split[2])
	except:
		error_exit_with_message("cannot convert FOLDER_split[2] to int, FOLDER_split[2]= %s, FOLDER= %s" %(FOLDER_split[2], FOLDER) )

	if(FOLDER in full_length_FOLDER_list): continue

	partial_length_FOLDER_list.append(FOLDER)
	

print
print "partial_length_FOLDER_list=", partial_length_FOLDER_list
print

PARTIAL_LENGTH_REGION_FOLDER="PARTIAL_LENGTH_REGION_FOLDER/"

if(exists(PARTIAL_LENGTH_REGION_FOLDER)): error_exit_with_message("PARTIAL_LENGTH_REGION_FOLDER (%s) already exist!" %(PARTIAL_LENGTH_REGION_FOLDER) )
submit_subprocess("mkdir %s" %(PARTIAL_LENGTH_REGION_FOLDER))

for partial_FOLDER in partial_length_FOLDER_list:
	start_location=partial_FOLDER
	final_location=PARTIAL_LENGTH_REGION_FOLDER + '/' + partial_FOLDER

	if(exists(start_location)==False):error_exit_with_message("start_location (%s) doesn't exist!" %(start_location) )
	if(exists(final_location)==True):error_exit_with_message("final_location (%s) already exist!" %(final_location) )


	submit_subprocess("mv %s %s " %(start_location, final_location))

print "JOB SUCESSFULLY COMPLETED"


