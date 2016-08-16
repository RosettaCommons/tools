#!/usr/bin/env python

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
from os import system
from os.path import basename, dirname, exists, expanduser, abspath
from time import sleep
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_parse_options import parse_options
######################################################################


def get_full_length_REGION_FOLDER_list(final_cluster_dag_file, verbose=False):

	if(exists(final_cluster_dag_file)==False): error_exit_with_message("final_cluster_dag_file (%s) doesn't exist! " %(final_cluster_dag_file) ) 

	if(verbose):
		print
		print "final_cluster_dag_file= ", final_cluster_dag_file
		print


	dag_string_list=open( final_cluster_dag_file).readlines()

	num_argument_line_found=0
	
	dag_args_list=[]

	for dag_string_args in dag_string_list:

		curr_dag_args_list=dag_string_args.split()

		if(curr_dag_args_list[0]=="arguments"): 
			dag_args_list=curr_dag_args_list
			num_argument_line_found+=1
			
			
	if(num_argument_line_found!=1): error_exit_with_message('num_argument_line_found!=1' )

	if(verbose):
		print
		print "dag_args_list= ", dag_args_list
		print
		


	full_length_silent_file_list= parse_options( dag_args_list, "in:file:silent", [""])


	if(len(full_length_silent_file_list[0])==0 and len(full_length_silent_file_list)==1): 
		error_exit_with_message('cannot find "in:file:silent"!' )

	if(verbose):
		print
		print "full_length_silent_file_list=", full_length_silent_file_list
		print



	full_length_FOLDER_list=[]

	for silent_file in full_length_silent_file_list:
		full_length_FOLDER=silent_file.replace("_sample.cluster.out","").upper()
		if(exists(full_length_FOLDER)==False):  error_exit_with_message("full_length_FOLDER (%s) doesn't exist!" %(full_length_FOLDER) )
		full_length_FOLDER_list.append(full_length_FOLDER)


	if(verbose):
		print
		print "full_length_FOLDER_list=", full_length_FOLDER_list
		print

	return full_length_FOLDER_list




