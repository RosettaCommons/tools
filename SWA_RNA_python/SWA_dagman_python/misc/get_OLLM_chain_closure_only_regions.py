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


from SWA_trace_pathway_util import *
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options


#get_OLLM_chain_closure_only_regions.py -num_elements 5


num_elements= parse_options( argv, "num_elements", 0)
folder_layer= parse_options( argv, "folder_layer", 1)

if(num_elements==0): error_exit_with_message("user need to pass in positive integer for num_elements!" )

All_region_file_list=[]

for j in range(num_elements):

	i=(j+1) % num_elements

	folder_name="REGION_%d_%d/" % (i,j)

	if(exists(folder_name)==False): error_exit_with_message("folder_name (%s) doesn't exist! " %(folder_name) ) 

	print "folder_name= %s" %(folder_name)

	region_file_list=[]

	if(i==0 or j==0):
		region_file_list = glob( "%s/start_from_region_*_*_sample_filtered.out" %(folder_name))	
	else:
		region_file_list = glob( "%s/start_from_region_*_*_and_*_*_sample_filtered.out" %(folder_name))

	All_region_file_list.extend(region_file_list)

for n in range(len(All_region_file_list)):
	region_file=All_region_file_list[n]

	for ii in range(folder_layer):
		region_file="../" + region_file

	All_region_file_list[n]=region_file


print "All_region_file_list= %s" %( list_to_string(All_region_file_list) )
