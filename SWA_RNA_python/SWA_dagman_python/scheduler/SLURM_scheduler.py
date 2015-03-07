#!/usr/bin/env python

######################################################################
from SWA_dagman_python.utility.SWA_util import *

import string
import subprocess
from glob import glob 

MPI_SCHEDULER=False

SCHEDULER_TYPE="SLURM_scheduler"

######################################################################
# Use this to run cluster using a SLURM (Simple Linux Utility for Resource Management) scheduler.
# Note: current setup requires some PBS emulation commands (i.e., 'qstat -f')

# Note this implementation is somewhat specific to the Stanford's SHERLOCK
# clusters (see below).
######################################################################

def get_job_id_list():
	job_id_list = []
	for FILE in ['JOB_ID.txt'] + glob('SLAVE_JOBS/*/JOB_ID.txt'):
		with open(FILE, 'r') as f:
			job_id_list += string.split(string.strip(f.read()),'\n')
	job_id_list = [string.strip(id) for id in job_id_list if len(string.strip(id))]
	return job_id_list

def get_job_ids_str(delimiter=' '):
	return string.join(get_job_id_list(), str(delimiter))

######################################################################

def job_status_command():
	command = "squeue"
	command += " --format='%i %j %u %t %B'"
	command += " --sort='i'"
	command += " --jobs=" + get_job_ids_str(delimiter=',') 
	return command

def kill_job_command():
	return 'scancel'  #NOTE: scancel --user=user_name.

def running_job_state_str():
	return 'R'

def pending_job_state_str():
	return 'PD' #job is queued, eligible to run or routed.
	#Note another possibility is 'S' #Job is suspended.

######################################################################
def get_temp_qstat_filename(): #Wrapper that call this function have the duty to delete temp_data_filename after using it!

	trail_num=0

	while(True):

		trail_num+=1

		if(trail_num==10): error_exit_with_message("Unsuccessfully request job_status after 10 trails")

		temp_data_filename="temp_job_stat_V%s.txt" %(str(trail_num).zfill(3) )

		#if(exists(temp_data_filename)): submit_subprocess_allow_retry("rm %s" %(temp_data_filename))

		if(exists(temp_data_filename)): continue

		submit_subprocess_allow_retry( job_status_command() + ' > %s' %(temp_data_filename))

		return temp_data_filename

######################################################################
##May 09, 2012 HACKY for debugging purposes
def update_scheduler_queue_log_command():

	temp_data_filename=get_temp_qstat_filename()

	update_scheduler_queue_log_filename="LOG_DEBUG_QSTAT_FULL.txt"

	print "updating scheduler_queue_log (%s)" %(update_scheduler_queue_log_filename)

	submit_subprocess_allow_retry("cp %s %s" %(temp_data_filename, update_scheduler_queue_log_filename))

	submit_subprocess_allow_retry("rm %s" %(temp_data_filename))


######################################################################
def kill_all_queued_jobs_with_prefix(job_prefix, verbose=True):

	if(verbose): print "kill_all_queued_jobs_with_prefix(prefix=%s)" %(job_prefix)

	sleep(10) #April 21, 2012: Give sometime for dead jobs to properly exit.

	if(isinstance(job_prefix, str)==False):
		print "ERROR: job_prefix=", job_prefix
		error_exit_with_message("job_prefix object is not a string!")


	jobs_stat_list = get_queued_jobs_status() #This gives the job_info_list instead of the string.

	num_job_killed=0

	for job_ID in range(len(jobs_stat_list)):

		job_info=copy.deepcopy(jobs_stat_list[job_ID])

		job_name=job_info['JOB_NAME']

		if(job_name[:len(job_prefix)]==job_prefix):
			num_job_killed+=1

			job_ID=job_info['JOBID']

			command=kill_job_command()
			command+=' %s' %(job_ID)

			if(verbose): print command
			submit_subprocess_allow_retry( command )

	if(num_job_killed==0):
		print "WARNING: num_job_killed=0 for job_prefix (%s)" %(job_prefix)
		print
		#print "ERROR: job_prefix=%s" %(job_prefix)
		#print "ERROR: num_job_with_matching_prefix=%d" %(num_job_with_matching_prefix)
		#error_exit_with_message("num_job_with_matching_prefix=0")

	return num_job_killed

######################################################################
def 	master_kill_all_slave_jobs_and_exit_scheduler_specific(exit_message):

	JOBDIR = 'SLAVE_JOBS/'
	job_prefix = abspath( JOBDIR ).replace('/','_')[1:]

	sys.stdout.flush()
	sys.stderr.flush()

 	kill_all_queued_jobs_with_prefix(job_prefix, verbose=True)

	sys.stdout.flush()
	sys.stderr.flush()
	exit_message2="Killed all slave jobs and about to exit master script, " + exit_message
	error_exit_with_message(exit_message2)


######################################################################

def get_queued_jobs_status():

	job_info_list = []
	temp_data_filename = get_temp_qstat_filename()
	data = safe_open(temp_data_filename, 'r')   #Sept 30, 2010

	'''
	SHERLOCK: 
	$ squeue --format="%i %j %u %t %B" --sort="i" --jobs=1713223,1713224,1713125
	JOBID NAME USER ST EXEC_HOST
	1713125 l2_viral_rna_pseudoknot_MASTER geniesse R sh-3-20
	1713223 l2_viral_rna_pseudoknot_SLAVE_JOBS_97 geniesse R sh-2-29
	1713224 l2_viral_rna_pseudoknot_SLAVE_JOBS_98 geniesse PD n/a
	'''

	# init empty header column list
        H = []
	
	# each line holds info for a single jobid
	for idx, line in enumerate(data.readlines()):
			
		# split into columns
		cols = string.split(string.strip(line))

		# skip empty lines, okay for now
		if not len(cols):
			continue
		
		# read header 
		if not len(H):
			H = cols
			continue
		
		# init new job_info
		job_info = {}

		# set values of attributes, based on index in header
		job_info['JOBID']     = cols[H.index('JOBID')]
		job_info['JOB_NAME']  = cols[H.index('NAME')]
		job_info['USER']      = cols[H.index('USER')]
		job_info['STATE']     = cols[H.index('ST')]
		job_info['EXEC_HOST'] = cols[H.index('EXEC_HOST')]

		# append job_info to job_info_list
		job_info_list.append(job_info)
		
	submit_subprocess_allow_retry("rm %s" %(temp_data_filename))

	return job_info_list


######################################################################
def queue_job_command(job_name, outfile, errfile, job_script, job_dir_name, walltime, memory_limit, memory_reserve, queue_name):

	#walltime (is hours): the job wall-clock run time limit.
	#memory_limit (in MB): job is killed when it exceeds this memory limit. This parameter does not guarantee memory allocation, it's just a threshold. (#NOT IMPLEMENTED for PBS scheduler)
	#memory_reserve (in MB): Amount of memory reserved/allocated for this job. This parameter guarantee memory allocation [Note memory_reserve should be less-than-or-equal-to memory_limit].
	#Assume using default queue. (For SHERLOCK cluster, this is the 'normal' queue | #SBATCH -p normal)

	### current walltime limit for SHERLOCK is 48 hours
	if walltime > 48: walltime = 48

	sbatch_submit_file="JOB_SCRIPT.sbatch" #"sbatch_submit_file"
	sbatch_outfile = outfile + "_SBATCH" 
	sbatch_errfile = errfile + "_SBATCH" 

	if(job_dir_name!=""): 
		sbatch_submit_file = job_dir_name + "/" + sbatch_submit_file
		if not '/' in sbatch_outfile:
			sbatch_outfile = job_dir_name + "/" + sbatch_outfile
		if not '/' in sbatch_outfile:
			sbatch_errfile = job_dir_name + "/" + sbatch_errfile

	if(exists(sbatch_submit_file)): submit_subprocess("rm %s" %(sbatch_submit_file))
	if(exists(sbatch_outfile)): submit_subprocess("rm %s" %(sbatch_outfile))
	if(exists(sbatch_errfile)): submit_subprocess("rm %s" %(sbatch_errfile))

	with open( sbatch_submit_file, 'w') as SBATCH_JOB_SCRIPT:

		SBATCH_JOB_SCRIPT.write( '#!/bin/bash\n\n' )

		SBATCH_JOB_SCRIPT.write( '#SBATCH -J %s\n' %(job_name))
		SBATCH_JOB_SCRIPT.write( '#SBATCH -o %s\n' %(sbatch_outfile))
		SBATCH_JOB_SCRIPT.write( '#SBATCH -e %s\n' %(sbatch_errfile))
		SBATCH_JOB_SCRIPT.write( '#SBATCH -t %d:00:00\n' %(walltime))
		SBATCH_JOB_SCRIPT.write( '#SBATCH -n 1\n' )
		SBATCH_JOB_SCRIPT.write( '#SBATCH -N 1\n' )

		if(memory_reserve!=0):
			SBATCH_JOB_SCRIPT.write( '#SBATCH --mem=%d\n' %(memory_reserve))  #memory per job
			SBATCH_JOB_SCRIPT.write( '#SBATCH --mem-per-cpu=%d\n' %(memory_reserve)) #memory per CPU process. Note mem and mem-per-cpu is the same for a single CPU job (which is the case here!).

		if(queue_name!="DEFAULT"): #Note on Biox2, defualt is the SP queue.
			SBATCH_JOB_SCRIPT.write( '#SBATCH -p %s\n' %(queue_name))

		
		current_directory=os.path.abspath(".")

		SBATCH_JOB_SCRIPT.write( '\ncd %s\n\n' %(current_directory)) #Important since by default SBATCH job starts in the home directory.

		#SBATCH_JOB_SCRIPT.write( 'cd $SLURM_SUBMIT_DIR\n\n' ) #This doesn't work if the node submitting the job is not the frontend (i.e. is itself a working_node)

		JOB_ID_FILE = 'JOB_ID.txt'
		
		if(job_dir_name!=""):
			JOB_ID_FILE = job_dir_name + '/' + JOB_ID_FILE

		SBATCH_JOB_SCRIPT.write( 'echo $SLURM_JOBID > %s\n' % JOB_ID_FILE )

		SBATCH_JOB_SCRIPT.write( '\n%s >%s 2>%s\n' %(job_script, outfile, errfile))
		
	submit_subprocess_allow_retry('sbatch %s' %(sbatch_submit_file))


######################################################################
def submit_DAG_job_scheduler_specific(master_wall_time, master_memory_reserve, num_slave_nodes, dagman_file, verbose):


	DAG_master_script_EXE=get_PYEXE("dagman/DAG_continuous.py")

	master_tag = abspath( "MASTER" ).replace("/","_")[1:]

	master_outfile='master_log.out'
	master_errfile='master_log.err'

	if( exists(master_outfile) ): submit_subprocess("rm %s" %(master_outfile))
	if( exists(master_errfile) ): submit_subprocess("rm %s" %(master_errfile))

	if( exists("JOB_ID.txt") ): submit_subprocess("rm -rf JOB_ID.txt")
	if( exists("SLAVE_JOBS/") ): submit_subprocess("rm -rf SLAVE_JOBS/")

	#Make slave_nodes wall_time 1 hour less than the master_node to ensure that the slave_nodes will die before the master_node!
	slave_wall_time=master_wall_time-1

	DAG_master_script_command=" %s -num_slave_jobs %d -slave_wall_time %d -dagman_file %s " %(DAG_master_script_EXE, num_slave_nodes, slave_wall_time, dagman_file)

	job_dir_name=""

	master_queue_name="normal" #This is specific to the SHERLOCK cluster (update if rhiju buys in).

	#Submit the master_job. #In this non-MPI implementation, the master_job will then submit the slave_jobs inside the DAG_continuous.py script.
	queue_job_command(master_tag, master_outfile, master_errfile, DAG_master_script_command, job_dir_name, master_wall_time, master_memory_reserve, master_memory_reserve, master_queue_name)












