#!/usr/bin/env python

######################################################################
from SWA_dagman_python.utility.SWA_util import *



MPI_SCHEDULER=False

SCHEDULER_TYPE="SLURM_scheduler"

######################################################################
# Use this to run cluster using a SLURM (Simple Linux Utility for Resource Management) scheduler.
# Note: current setup requires some PBS emulation commands (i.e., 'qstat -f')

# Note this implementation is somewhat specific to the Stanford's SHERLOCK
# clusters (see below).
######################################################################
def queue_status_user_command():
	from os.path import expandvars
	return 'squeue --user=%s' % expandvars('$USER') 

def get_job_id_list():
	'''
	 		 JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
           1710723       dev     sdev geniesse  R      31:02      1 sh-5-35
    '''
	import string
	from subprocess import Popen, PIPE
	cmd = string.split(queue_status_user_command())
	out, err = Popen(cmd, stdout=PIPE, stderr=PIPE).communicate()
	for line in out.split('\n'):
		cols = line.split()
		if len(cols) < 1: continue
		if 'JOBID' in cols[0]: continue
		try:
			is_int = int(cols[0])
		except:
			continue
		job_id_list.append(str(cols[0]))
	return job_id_list

######################################################################

def job_status_command():
	return 'qstat -f %s' % (' '.join(get_job_id_list()))

def kill_job_command():
	return 'scancel'  #NOTE: scancel --user=user_name.

def running_job_state_str():
	return 'R'

def pending_job_state_str():
	return 'Q' #job is queued, eligible to run or routed.

	#Note another possibility is 'H' #Job is held.


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

	if(isinstance(job_prefix, str )==False):
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

	'''qstat -f:
	
	SHERLOCK:
	Job Id:	1710705
		Job_Name = STAR_Mapping%j
		Job_Owner = klane@sh-2-34
		job_state = Q
		queue = mcovert
		qtime = Fri Mar  6 15:18:08 2015
		mtime = Sat Mar  7 15:22:10 2015
		Account_Name = mcovert
		Priority = 1073
		euser = klane(6753)
		egroup = mcovert(38308)
		Resource_List.walltime = 24:00:00
		Resource_List.nodect = 1
		Resource_List.ncpus = 1

	BIOX3
	Job Id: 2323754.biox3-frontend-1.stanford.edu
	    Job_Name = biox3_home_geniesse_src_stepwise_benchmark_new_swa_revival_swa_
		ss_loop_mode_vdw_rep_screen_tether_jump_false_23s_rrna_2003_2012_SLAVE
		_JOBS_499
	    Job_Owner = geniesse@biox3-1-1.Stanford.EDU
	    job_state = Q
	    queue = SP
	    server = biox3-frontend-1.stanford.edu
	    Checkpoint = u
	    ctime = Fri Mar  6 09:00:57 2015
	    Error_Path = biox3-1-1.stanford.edu:/biox3/home/geniesse/src/stepwise_benc
		hmark/new/swa_revival/swa_ss_loop_mode_vdw_rep_screen_tether_jump_fals
		e/23s_rrna_2003_2012/SLAVE_JOBS/499/slave_jobs.err_QSUB
	    Hold_Types = n
	    Join_Path = n
	    Keep_Files = n
	    Mail_Points = a
	    mtime = Fri Mar  6 09:00:57 2015
	    Output_Path = biox3-1-1.stanford.edu:/biox3/home/geniesse/src/stepwise_ben
		chmark/new/swa_revival/swa_ss_loop_mode_vdw_rep_screen_tether_jump_fal
		se/23s_rrna_2003_2012/SLAVE_JOBS/499/slave_jobs.out_QSUB
	    Priority = 0
	    qtime = Fri Mar  6 09:00:57 2015
	    Rerunable = True
	    Resource_List.mem = 2048mb
	    Resource_List.nodect = 1
	    Resource_List.nodes = 1:ppn=1
	    Resource_List.pmem = 2048mb
	    Resource_List.walltime = 71:00:00
	    Variable_List = PBS_O_QUEUE=SP,PBS_O_HOME=/home/geniesse,

	'''


	temp_data_filename=get_temp_qstat_filename()

	data=safe_open(temp_data_filename, 'r')   #Sept 30, 2010

	#	JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
	job_info_list=[]

	while(True):

		line=data.readline()

		if(line==''): break #End of file!

		if(line[0:7]=="Job Id:"):

			JOB_OWNER='BLAH_BLAH_BLAH'
			JOB_NAME='BLAH_BLAH_BLAH'

			job_info={}
			job_info['JOBID']=line.split()[-1]

			line=data.readline()

			####################################################################################################
			if(line.split()[0]=='Job_Name'): #Parse Job_Name then Job_Owner
				JOB_NAME=line.split()[-1]

				while(True):
					line=data.readline()

					#print "DEBUG line.split()=", line.split()

					if(len(line.split())==0): line=data.readline() #April 22, 2012: Deal with special empty line case (see right below)

					if(len(line.split())!=1): break

					JOB_NAME+=line.split()[0]

					'''
					Job_Name = home_sripakpa_minirosetta_April_22_4X_CHEM_SHIFT_minimize_bench
					mark_run_main_folder_FARFAR_denovo_4X4_R2_Retrotransposon_SLAVE_JOBS_0
		      											<--Empty line here!
		  			Job_Owner = sripakpa@node-9-34.local
					'''

				if(line.split()[0]!='Job_Owner'): error_exit_with_message("line.split()[0]!='Job_Owner' for line=%s" %(line))
				if(line.split()[1]!='='): error_exit_with_message("line.split()[1]!='=' for line=%s" %(line))

				JOB_OWNER=line.split()[-1]


			elif(line.split()[0]=='Job_Owner'):  #Parse Job_Name then Job_Owner

				JOB_OWNER=line.split()[-1]

				line=data.readline()

				if(line.split()[0]!='Job_Name'): error_exit_with_message("line[0:9]!='Job_Owner' for line=%s" %(line))

				JOB_NAME=line.split()[-1]

				while(True):
					line=data.readline()

					#print "DEBUG line.split()=", line.split()

					if(len(line.split())==0): line=data.readline() #April 22, 2012: Deal with special empty line case.

					if(len(line.split())!=1): break

					JOB_NAME+=line.split()[0]

				if(line.split()[1]!='='): error_exit_with_message("line.split()[1]!='=' for line=%s" %(line))

			else:
				error_exit_with_message("Invalid first job status line=%s" %(line))
			####################################################################################################

			job_info['JOB_NAME']=JOB_NAME
			job_info['USER']=JOB_OWNER


			while(True):
				line=data.readline()

				if(line.split()[0]=='job_state'):
					job_info['STATE']=line.split()[-1]
					break


			job_info['EXEC_HOST']="NONE"

			if(job_info['STATE']==running_job_state_str()):

				check_line_list=[]

				while(True):
					line=data.readline()
					check_line_list.append(line)
					if(line=="") or (len(line)==0): error_exit_with_message("Reach end of file: incomplete parsing of qstat!")

					if(len(check_line_list)>30):
						for n in range(len(check_line_list)):
							print "ERROR #n: %s" %(check_line_list[n]),
						error_exit_with_message("len(check_line_list)>30!")

					#if(line.split('=')[0]=='exec_host'):
					if ( 'exec_host' in line ):
						job_info['EXEC_HOST']=line.split()[-1]
						break

			job_info_list.append(job_info)

	data.close()

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

	sbatch_submit_file="JOB.sbatch" #"sbatch_submit_file"

	if(job_dir_name!=""): sbatch_submit_file= job_dir_name + "/" + sbatch_submit_file

	if(exists(sbatch_submit_file)): submit_subprocess("rm %s" %(sbatch_submit_file))

	SBATCH_JOB = open( sbatch_submit_file, 'w')

	SBATCH_JOB.write( '#!/bin/bash\n\n' )

	SBATCH_JOB.write( '#SBATCH -J %s\n\n' %(job_name))
	SBATCH_JOB.write( '#SBATCH -o %s\n\n' %(outfile + "_SBATCH" ))
	SBATCH_JOB.write( '#SBATCH -e %s\n\n' %(errfile + "_SBATCH" ))
	SBATCH_JOB.write( '#SBATCH -t %d:00:00\n\n' %(walltime))
	SBATCH_JOB.write( '#SBATCH -n 1\n\n' )
	SBATCH_JOB.write( '#SBATCH -N 1\n\n' )


	if(memory_reserve!=0):
		SBATCH_JOB.write( '#SBATCH --mem=%d\n\n' %(memory_reserve))  #memory per job
		SBATCH_JOB.write( '#SBATCH --mem-per-cpu=%d\n\n' %(memory_reserve)) #memory per CPU process. Note mem and mem-per-cpu is the same for a single CPU job (which is the case here!).

	if(queue_name!="DEFAULT"): #Note on Biox2, defualt is the SP queue.
		SBATCH_JOB.write( '#SBATCH -p %s\n\n' %(queue_name))

	
	current_directory=os.path.abspath(".")

	SBATCH_JOB.write( 'cd %s\n\n' %(current_directory)) #Important since by default SBATCH job starts in the home directory.

	#SBATCH_JOB.write( 'cd $SLURM_SUBMIT_DIR\n\n' ) #This doesn't work if the node submitting the job is not the frontend (i.e. is itself a working_node)

	if(job_dir_name!=""):
		SBATCH_JOB.write( 'echo $SLURM_JOBID > %s/JOB_ID.txt\n\n' %(job_dir_name))
	else:
		SBATCH_JOB.write( 'echo $SLURM_JOBID > JOB_ID.txt\n\n')

	SBATCH_JOB.write( '%s >%s 2>%s\n\n' %(job_script, outfile, errfile))

	SBATCH_JOB.close()

	submit_subprocess_allow_retry('sbatch %s' %(sbatch_submit_file))

######################################################################
def submit_DAG_job_scheduler_specific(master_wall_time, master_memory_reserve, num_slave_nodes, dagman_file, verbose):


	DAG_master_script_EXE=get_PYEXE("dagman/DAG_continuous.py")

	master_tag = abspath( "MASTER" ).replace("/","_")[1:]

	master_outfile='master_log.out'
	master_errfile='master_log.err'

	if( exists(master_outfile) ): submit_subprocess("rm %s" %(master_outfile))
	if( exists(master_errfile) ): submit_subprocess("rm %s" %(master_errfile))

	if( exists("SLAVE_JOBS/") ): submit_subprocess("rm -rf SLAVE_JOBS/")

	#Make slave_nodes wall_time 1 hour less than the master_node to ensure that the slave_nodes will die before the master_node!
	slave_wall_time=master_wall_time-1

	DAG_master_script_command=" %s -num_slave_jobs %d -slave_wall_time %d -dagman_file %s " %(DAG_master_script_EXE, num_slave_nodes, slave_wall_time, dagman_file)

	job_dir_name=""

	master_queue_name="normal" #This is specific to the SHERLOCK cluster (update if rhiju buys in).



	#Submit the master_job. #In this non-MPI implementation, the master_job will then submit the slave_jobs inside the DAG_continuous.py script.
	queue_job_command(master_tag, master_outfile, master_errfile, DAG_master_script_command, job_dir_name, master_wall_time, master_memory_reserve, master_memory_reserve, master_queue_name)












