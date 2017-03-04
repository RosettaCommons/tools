#!/usr/bin/env python

######################################################################
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.utility.DAGMAN_util import create_generic_README_SUB

######################################################################

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
import os
from os.path import basename, dirname, exists, expanduser, abspath
from time import sleep
######################################################################

README_SUB_PY="README_SUB.py"
README_SETUP_PY="README_SETUP.py"
######################################################################

def get_num_job_queued_in_account():

	###Ensure the slave jobs is really dead!###
	bjobs_lines=popen_and_readlines("bjobs -w | sort -nk1 ", Is_master=False, tag="get_bjob_status.txt")

	return ( len(bjobs_lines) -1 )
######################################################################


def get_num_existing_slave_jobs():

	total_slave=0
	num_running=0
	num_pending=0

	slave_tag = abspath( "SLAVE_JOBS/" ).replace('/','_')
	#print "slave_tag=%s" %(slave_tag)

	###Ensure the slave jobs is really dead!###
	bjobs_lines=popen_and_readlines("bjobs -w | sort -nk1 ", Is_master=False, tag="get_bjob_status.txt")

	for line in bjobs_lines:
		
		if(line.find( slave_tag ) > 0): 
			total_slave+=1

			if(line.split()[2]=="RUN"): 
				num_running+=1
			elif(line.split()[2]=="PEND"): 
				num_pending+=1
			elif(line.split()[2]=="STAT"):
				pass
			else:
				error_exit_with_message("Invalid status for job=%s" %(line))		
							
	return (total_slave, num_running, num_pending)

######################################################################
def get_master_job_status():

	master_tag = abspath( "MASTER" ).replace("/","_")
	#print "master_tag=%s" %(master_tag)

	###Ensure the slave jobs is really dead!###
	bjobs_lines=popen_and_readlines("bjobs -w | sort -nk1 ", Is_master=False, tag="get_bjob_status.txt")

	master_job_status="NOT_QUEUED"

	num_master_job_line_found=0

	for line in bjobs_lines:
		
		if(line.find( master_tag ) > 0): 
			num_master_job_line_found+=1

			if(line.split()[2]=="RUN"): 
				master_job_status="RUNNING"
			elif(line.split()[2]=="PEND"): 
				master_job_status="PENDING"
			elif(line.split()[2]=="STAT"):
				pass
			else:
				error_exit_with_message("Invalid status for job=%s" %(line))		

	if(num_master_job_line_found>1): error_exit_with_message("num_master_job_line_found>1")
	
	return master_job_status
######################################################################

def check_job_is_done():

	if(exists("master_log.out")==False): return False

	master_log_out_tails=popen_and_readlines("less master_log.out | tail -n50", Is_master=False, tag="master_log_out_tails.txt")

	found_master_job_done_line=False

	for line in master_log_out_tails:
		if(line=="Successfuly RAN DAG_continuous.py!!\n"): found_master_job_done_line=True

	return found_master_job_done_line


######################################################################

def submit_job():

	print "python %s " %(README_SUB_PY)
	print "python %s " %(README_SETUP_PY)
	if(exists(README_SUB_PY)==False): error_exit_with_message("README_SUB_PY (%s) doesn't exist!" %(README_SUB_PY) )
	if(exists(README_SETUP_PY)==False): error_exit_with_message("README_SETUP_PY (%s) doesn't exist!" %(README_SETUP_PY))
	submit_subprocess("python %s " %(README_SUB_PY))
	submit_subprocess("python %s " %(README_SETUP_PY))

######################################################################

def attempt_resubmit_job():

	if(exists("master_log.err")==False): return False

	if(exists("master_log.out")==False): error_exit_with_message("master_log.err exist not master_log.out does not!")

	master_log_err_lines=safe_readlines("master_log.err")

	job_resubmitted=False

	if(len(master_log_err_lines)>0): #OK JOB died due to error. Try REsubmitting!
		print "NON-empty master_log.err!"
		for n in range(len(master_log_err_lines)):
			print "master_log.err line #%d: %s" %(n, master_log_err_lines[n]),

		#if(exists("SLAVE_JOBS/")==False): error_exit_with_message("SLAVE_JOBS/ doesn't exist!")

		print "---REsubmitting job!---" 

		#########################################################
		count=0
		while(True):
			count+=1
			OLD_folder="AUTO_OLD_%d" %(count)

			if(exists(OLD_folder)==False):
				print "moving OLD log files to %s" %(OLD_folder)
				submit_subprocess("mkdir %s" %(OLD_folder))

				if(exists("master_log.out")):	submit_subprocess("mv master_log.out %s/master_log.out" %(OLD_folder))
				if(exists("master_log.err")):	submit_subprocess("mv master_log.err %s/master_log.err" %(OLD_folder))
				if(exists("SLAVE_JOBS/")): 		submit_subprocess("mv SLAVE_JOBS/ %s/SLAVE_JOBS/" %(OLD_folder))
				if(exists("CONDER/")):				submit_subprocess("mv CONDER/ %s/CONDER/" %(OLD_folder))
				if(exists("COMMON_ARGS/")):		submit_subprocess("mv COMMON_ARGS/ %s/COMMON_ARGS/" %(OLD_folder))
				if(exists("KEEP_LOG_FILE/")):	submit_subprocess("cp KEEP_LOG_FILE %s/KEEP_LOG_FILE" %(OLD_folder))
				break

		#########################################################
		sleep( 2 ) #SLAVE jobs that haven't died yet should start dying after moving the SLAVE_JOBS folder.

		(total_slave,  num_running_slave, num_pending_slave)=get_num_existing_slave_jobs()			

		master_node_status = get_master_job_status()

		if(total_slave>0): error_exit_with_message("NON-empty master_log.err BUT total_slave_nodes (%s) STILL EXIST in queue!" %(total_slave) ) 

		if(master_node_status!="NOT_QUEUED"): error_exit_with_message("NON-empty master_log.err BUT master_job_is_running is STILL QUEUED!") 

		#########################################################
		job_resubmitted=True
		submit_job()

	return job_resubmitted

######################################################################


def attempt_submit_job():

	master_node_status=get_master_job_status()

	job_submitted=False

	if(master_node_status=="NOT_QUEUED"): #JOB have not been previosly submitted (SO DOES NOT EXIST IN QUEUE)

		if(exists("master_log.out")): error_exit_with_message("master_log.out already exist!")
		if(exists("master_log.err")): error_exit_with_message("master_log.err already exist!")
		if(exists("SLAVE_JOBS/")): error_exit_with_message("SLAVE_JOBS/ already exist!")

		job_submitted=True
		submit_job()

	elif(master_node_status=="PENDING"): #JOB have already been previously submitted BUT PENDING

		if(exists("master_log.out")): error_exit_with_message("master_log.out already exist!")
		if(exists("master_log.err")): error_exit_with_message("master_log.err already exist!")
		if(exists("SLAVE_JOBS/")): error_exit_with_message("SLAVE_JOBS/ already exist!")

	elif(master_node_status=="RUNNING"):

		(total_slave,  num_running_slave, num_pending_slave)=get_num_existing_slave_jobs()		
		#if(total_slave==0): error_exit_with_message("Job should be running but total_slave==0!") 
		print "job is already RUNNING | total_slave=%4s | num_running=%4s | num_pending=%4s " %(total_slave, num_running_slave, num_pending_slave)	
	##################################################################################################################
	else:
		error_exit_with_message("Invalid master_node_status (%s)" %(master_node_status))

	return job_submitted

###########################################################################

def check_status():

	(total_slave, num_running_slave, num_pending_slave)=get_num_existing_slave_jobs()

	master_node_status=get_master_job_status()

	job_is_done=check_job_is_done()

	if(exists("master_log.err")): #Check for possible errors.
		master_log_err_lines=safe_readlines("master_log.err")

		if(len(master_log_err_lines)>0):
			print "NON-empty master_log.err!"
			for n in range(len(master_log_err_lines)):
				print "master_log.err line #%d: %s" %(n, master_log_err_lines[n]),
			print "job DIED!      | total_slave=%4s | num_running=%4s | num_pending=%4s " %(total_slave, num_running_slave, num_pending_slave)			
			error_exit_with_message("	NON-empty master_log.err!")			


	if(job_is_done):		
		print "master_job is DONE!  | total_slave=%4s | num_running=%4s | num_pending=%4s " %(total_slave, num_running_slave, num_pending_slave)			
		if(total_slave>0): error_exit_with_message("job_is_done but total_slave_nodes (%s) STILL EXIST in queue!" %(total_slave) )
		return

	if(master_node_status=="NOT_QUEUED"):
		print "master_job IS NOT YET SUBMITTED"
	
	elif(master_node_status=="PENDING"):
		print "master_job IS PENDING "

	elif(master_node_status=="RUNNING"):

		if(exists("master_log.out")==False): 
			print "master_job is RUNNING but master_log.out does not exist!"
			return

		if(exists("master_log.err")==False): 
			print "master_job is RUNNING but master_log.err does not exist!"
			return

		print "master_job is RUNNING | total_slave=%4s | num_running=%4s | num_pending=%4s " %(total_slave, num_running_slave, num_pending_slave)			

	else:
		error_exit_with_message("Invalid master_node_status (%s)" %(master_node_status))

###########################################################################

