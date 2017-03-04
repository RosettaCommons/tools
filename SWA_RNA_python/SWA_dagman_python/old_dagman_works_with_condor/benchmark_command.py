#!/usr/bin/env python

######################################################################
from benchmark_command_util import *

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
from SWA_dagman_python.parser.SWA_parse_benchmark import parse_benchmark_job_file
######################################################################

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
import os
from os.path import basename, dirname, exists, expanduser
from time import sleep
######################################################################

#~/SWA_RNA_python/SWA_dagman_python/dagman/benchmark_command.py  -algorithm check_status -job_folder_patterns FARFAR_denovo/^/ SWA_denovo/^/ 


#~/SWA_RNA_python/SWA_dagman_python/dagman/benchmark_command.py -algorithm submit_jobs -job_folder_patterns FARFAR_denovo/^/ SWA_denovo/^/  -num_slave_nodes_per_job 100 >> log_submit_job_benchmark_command.out 2> log_submit_job_benchmark_command.err

#  


START_argv=copy.deepcopy(argv)

print "---------------------------------------------------------------------------------------------------------------------"
print "---------------------------------------------------------------------------------------------------------------------"
print "Enter %s" %(list_to_string(START_argv)) 
print

job_folder_patterns= parse_options( argv, "job_folder_patterns", [""] )

if(job_folder_patterns==""): error_exit_with_message('User need to specify job_folder_patterns!')


algorithm= parse_options( argv, "algorithm", "" )

if(algorithm==""): error_exit_with_message('User need to specify algorithm!')

if(algorithm not in ["submit_jobs", "REsubmit_jobs", "check_status"] ):  error_exit_with_message("Invalid algorithm: (%s)" %(algorithm)) 

num_slave_nodes_per_job=parse_options( argv, "num_slave_nodes_per_job", 100 )

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )

############################################################################################

main_folder=os.path.abspath("main_folder/")
if(exists(main_folder)==False): error_exit_with_message("main_folder (%s) doesn't exist!" %(main_folder)) 
os.chdir( main_folder)

job_folder_list=[]

print "job_folder_patterns=", job_folder_patterns

for pattern_ID in range(len(job_folder_patterns)):

	job_folder_patterns[pattern_ID]=job_folder_patterns[pattern_ID].replace("^", "*")

	job_folder_list.extend( glob( job_folder_patterns[pattern_ID] )	)

print "job_folder_list=", job_folder_list

print "------------------------------------------------------------------------------------------"
print "------------------------------------------------------------------------------------------"

num_job_queued_so_FAR=get_num_job_queued_in_account()



submitted_job_list=[]

num_rounds=1

if(algorithm=="submit_jobs"): 
	print "BEFORE submitting jobs, num_job_queued=%s" %(num_job_queued_so_FAR)
	num_rounds=2
else:
	print "num_job_queued=%s" %(num_job_queued_so_FAR)

for round_ID in range(num_rounds):

	if(algorithm=="submit_jobs" and round_ID==1): sleep(5)

	for job_folder in job_folder_list:

		print "---------------------------------------------------------------------------------------------------------"
		print "round_ID=%s | Current job_folder: " %(round_ID), job_folder 

		os.chdir( main_folder)
		os.chdir( job_folder )

		if(algorithm=="submit_jobs"):

			job_is_done=check_job_is_done()

			if(job_is_done): 
				if(round_ID==1): print "job is already done!"
				continue

			if(job_folder in submitted_job_list):
				print "job was just submitted in prior round_ID!"
				continue

			job_resubmitted=False

			#Check if error occur, if so that make sure that no leftover jobs [including master] then resubmit.
			if(round_ID==0): job_resubmitted=attempt_resubmit_job()

			#OK if no error, then check whether the job is already running, if not then submit it!
			if(round_ID==1): job_resubmitted=attempt_submit_job()


			if(job_resubmitted): 
				submitted_job_list.append(job_folder)
				num_job_queued_so_FAR+=num_slave_nodes_per_job

			if(num_job_queued_so_FAR>600):
				print "num_job_running_so_FAR=(%s)>600 EARLY EXIT!" %(num_job_queued_so_FAR)
				print "Submitted JOBS:",  submitted_job_list
				sys.exit(1)

		elif(algorithm=="check_status"):
	
			check_status()
			
		else: #FARFAR_easy_cat
			error_exit_with_message("Invalid algorithm=%s" %(algorithm))



print "------------------------------------------------------------------------------------------"
print "------------------------------------------------------------------------------------------"
print "SUCCESSFULLY RAN benchmark_command.py!"
