#!/usr/bin/env python

######################################################################
from SWA_dagman_python.utility.SWA_util import *

######################################################################
#Use this to run cluster using a LSF (Load Sharing Facility) scheduler
######################################################################

def job_status_command():
	return 'bjobs -w'

def kill_job_command():
	return 'bkill'

def job_state_column_name():
	return 'STAT'

def job_name_column_name():
	return 'JOB_NAME'

def job_ID_column_name():
	return 'JOBID'

def running_job_state_str():
	return 'RUN'

def pending_job_state_str():
	return 'PEND'

######################################################################
##May 09, 2012 HACKY for debugging PBS_scheduler !
def update_scheduler_queue_log_command():

	error_exit_with_message("This function is not yet implemented!")


def get_queued_jobs_status():

	jobs_stat_lines = popen_and_readlines(job_status_command(), tag="temp_job_stat.txt" )

	'''
	JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
	942172  kwipapat RUN   SP         node-1-16.local node-1-19   _home_kwipapat_minirosetta_Feb_06_Benchmark_hairpin_loop_rebuild_FULL_LENGTH_AND_FINAL_BULGE_main_folder_CUUG_tetraloop_SLAVE_JOBS_80 Feb  8 00:52
	942148  kwipapat RUN   SP         node-1-16.local node-2-26   _home_kwipapat_minirosetta_Feb_06_Benchmark_hairpin_loop_rebuild_FULL_LENGTH_AND_FINAL_BULGE_main_folder_CUUG_tetraloop_SLAVE_JOBS_56 Feb  8 00:52
	942169  kwipapat RUN   SP         node-1-16.local node-6-14   _home_kwipapat_minirosetta_Feb_06_Benchmark_hairpin_loop_rebuild_FULL_LENGTH_AND_FINAL_BULGE_main_folder_CUUG_tetraloop_SLAVE_JOBS_77 Feb  8 00:52
	942848  kwipapat RUN   SP         node-4-15.local node-6-14   _home_kwipapat_minirosetta_Feb_06_Benchmark_hairpin_loop_rebuild_FULL_LENGTH_AND_FINAL_BULGE_main_folder_TRNA_YE_MET_SLAVE_JOBS_50 Feb  8 02:34
	942998  kwipapat RUN   SP         node-4-15.local node-6-14   _home_kwipapat_minirosetta_Feb_06_Benchmark_hairpin_loop_rebuild_FULL_LENGTH_AND_FINAL_BULGE_main_folder_TRNA_YE_MET_SLAVE_JOBS_200 Feb  8 02:34
	942140  kwipapat RUN   SP         node-1-16.local node-2-21   _home_kwipapat_minirosetta_Feb_06_Benchmark_hairpin_loop_rebuild_FULL_LENGTH_AND_FINAL_BULGE_main_folder_CUUG_tetraloop_SLAVE_JOBS_48 Feb  8 00:52
	942143  kwipapat RUN   SP         node-1-16.local node-8-8    _home_kwipapat_minirosetta_Feb_06_Benchmark_hairpin_loop_rebuild_FULL_LENGTH_AND_FINAL_BULGE_main_folder_CUUG_tetraloop_SLAVE_JOBS_51 Feb  8 00:52
	942165  kwipapat RUN   SP         node-1-16.local node-5-16   _home_kwipapat_minirosetta_Feb_06_Benchmark_hairpin_loop_rebuild_FULL_LENGTH_AND_FINAL_BULGE_main_folder_CUUG_tetraloop_SLAVE_JOBS_73 Feb  8 00:52
	'''

	return jobs_stat_lines

######################################################################
def queue_job_command(job_name, outfile, errfile, job_script, job_dir_name, walltime, memory_limit, memory_reserve):

	#walltime (is hours): the job wall-clock run time limit.
	#memory_limit (in MB): job is killed when it exceeds this memory limit. This parameter does not guarantee memory allocation, it's just a threshold.
	#memory_reserve (in MB): Amount of memory reserved/allocated for this job. This parameter guarantee memory allocation [Note memory_reserve should be less-than-or-equal-to memory_limit].
	#Assume using default queue. (For BIOX2 cluster, this is the SP qeuue)

	command='bsub -J %s -o %s -e %s -W %d:0 -M %d ' %(job_name, outfile, errfile, wall_time, memory_limit) 

	if(memory_reserve!=0): command += '-R "rusage[mem=%d] ' %(memory_reserve)

	print command
	submit_subprocess( command )

	'''
	if(SEPERATE_CLUSTERER_SLAVE_NODES and (n<NUM_CLUSTERER_SLAVE_NODE ) ): #HACKY. High memory node for clustering!
		queue_job(job_name, outfile, errfile, dir_name, job_script, walltime=168, memory_soft=0, memory_strict=0, verbose=True)

		command = 'bsub -W 140:0 -M 8192000 -R "rusage[mem=8192]" -o %s -e %s -J %s %s -slave_dir %s ' % (slave_outfile, slave_errfile, slave_tag, SLAVE_EXE, SLAVE_DIR)
		#command = 'bsub -W 140:0 -M 5000000 -R "rusage[mem=5000]" -o %s -e %s -J %s %s -slave_dir %s ' % (slave_outfile, slave_errfile, slave_tag, SLAVE_EXE, SLAVE_DIR)
		#command = 'bsub -W 140:0 -M 4000000 -o %s -e %s -J %s %s -slave_dir %s ' % (slave_outfile, slave_errfile, slave_tag , SLAVE_EXE, SLAVE_DIR)
	else: #Standard node.
		command = 'bsub -W 140:0 -M 4000000 -o %s -e %s -J %s %s -slave_dir %s ' % (slave_outfile, slave_errfile, slave_tag , SLAVE_EXE, SLAVE_DIR)
		#command = 'bsub -W 140:0 -M 5000000 -R "rusage[mem=5000]" -o %s -e %s -J %s %s -slave_dir %s ' % (slave_outfile, slave_errfile, slave_tag, SLAVE_EXE, SLAVE_DIR)
	'''

