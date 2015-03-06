#!/usr/bin/env python

######################################################################


from SWA_dagman_python.utility.SWA_util import *

######################################################################

import os

COMPUTER_CLUSTER_NAME=str(os.environ['COMPUTER_CLUSTER_NAME'])

print "ENVIROMENT_VARIABLE: COMPUTER_CLUSTER_NAME=%s" %(COMPUTER_CLUSTER_NAME)

if COMPUTER_CLUSTER_NAME == "LOCAL_TEST":
	# from LSF_scheduler import *
	from PBS_scheduler import *
	# from MPI_SGE_scheduler import *

elif COMPUTER_CLUSTER_NAME in ["BIOX2-STANFORD",
                               "BIOX3-STANFORD"]:

	from PBS_scheduler import *

elif COMPUTER_CLUSTER_NAME in ["LONESTAR-TACC-XSEDE",
                               "RANGER-TACC-XSEDE",
                               "STAMPEDE-TACC-XSEDE"]:


	from MPI_SGE_scheduler import *

elif COMPUTER_CLUSTER_NAME in ["TRESTLES-SDSC-XSEDE", 
							   "KRAKEN-NICS-XSEDE"]:

	from MPI_PBS_scheduler import *

elif COMPUTER_CLUSTER_NAME in ["SHERLOCK-STANFORD"]:

	from SLURM_scheduler import *

else:
	error_exit_with_message("Unsupported COMPUTER_CLUSTER_NAME=%s" %(COMPUTER_CLUSTER_NAME))



###Functions in scheduler_util are cluster/environment independent.
###Cluster/environment dependent command should be put in the specific scheduler python script (e.g. PBS_scheduler).

############################################################################
def update_scheduler_queue_log():

	update_scheduler_queue_log_command()

############################################################################
def queue_job(job_name, outfile, errfile, job_script, job_dir_name, walltime=24, memory_limit=4096, memory_reserve=0, queue_name="DEFAULT", verbose=True): 

	#job_dir_name is relative path to the python script's current directory.

	print "queue_job_command: job_name=%s | outfile=%s | errfile=%s | job_script=%s | job_dir_name=%s | walltime=%s | memory_limit=%s | memory_reserve=%s | queue_name =%s" %(job_name, outfile, errfile, job_script, job_dir_name, walltime, memory_limit, memory_reserve, queue_name)


	############################################################################
	if(isinstance( job_name, str )==False): 
		print "ERROR: job_name:", job_name
		error_exit_with_message("job_name object is not an str!")

	if(isinstance( outfile, str )==False): 
		print "ERROR: outfile:", outfile
		error_exit_with_message("outfile object is not an str!")

	if(isinstance( errfile, str )==False): 
		print "ERROR: errfile:", errfile
		error_exit_with_message("errfile object is not an str!")

	if(isinstance( job_script, str )==False): 
		print "ERROR: job_script:", job_script
		error_exit_with_message("job_script object is not an str!")

	if(isinstance( job_dir_name, str )==False): 
		print "ERROR: job_dir_name:", job_dir_name
		error_exit_with_message("job_dir_name object is not an str!")

	if(isinstance( walltime, int )==False): 
		print "ERROR: walltime:", walltime
		error_exit_with_message("walltime object is not an int!")

	if(isinstance( memory_limit, int )==False): 
		print "ERROR: memory_limit:", memory_limit
		error_exit_with_message("memory_limit object is not an int!")

	if(isinstance( memory_reserve, int )==False): 
		print "ERROR: memory_reserve:", memory_reserve
		error_exit_with_message("memory_reserve object is not an int!")

	if(isinstance( queue_name, str )==False): 
		print "ERROR: queue_name:", queue_name
		error_exit_with_message("queue_name object is not an str!")
	############################################################################


	start_dir=os.path.abspath(".")

	#if(job_dir_name!=""):
	#	if(exists(job_dir_name)==False): submit_subprocess("mkdir -p %s" %(job_dir_name))
	#	os.chdir(job_dir_name)

	print "start_dir=%s | job_dir_name=%s " %(start_dir, job_dir_name)

	queue_job_command(job_name, outfile, errfile, job_script, job_dir_name, walltime, memory_limit, memory_reserve, queue_name)

	os.chdir(start_dir)

####################################################################

def get_job_col_index(column_list, col_name):

	if(column_list.count(col_name)!=1): 
		print "ERROR: job_status column_list=" %(list_to_string(column_list)[1:])
		error_exit_with_message("column_list.count('%s')!=1" %(col_name))

	col_index=column_list.index(col_name)

	return col_index

######################################################################

def get_job_stat_colname_line():

	return "%s   %s   %s   %20s   %s" %('JOBID', 'JOB_NAME', 'USER', 'EXEC_HOST', 'STATE')

######################################################################

def job_info_to_job_stat_line(job_info):

	return "%s   %s   %s   %20s   %s" %(job_info['JOBID'], job_info['JOB_NAME'], job_info['USER'], job_info['EXEC_HOST'], job_info['STATE']	)


######################################################################

def get_queued_jobs_status_lines():	

	job_info_list=get_queued_jobs_status()

	job_stat_lines=[]

	job_stat_lines.append(get_job_stat_colname_line())

	for job_info in job_info_list:
		job_stat_lines.append(job_info_to_job_stat_line(job_info))

	return job_stat_lines


######################################################################

def print_queued_jobs_status():

	jobs_stat_lines = get_queued_jobs_status_lines()

	print "---------------------------------------queued_jobs_status:-------------------------------------------"
	for n in range(len(jobs_stat_lines)):
		print "#%4d: %s" %(n, jobs_stat_lines[n])
	print "-----------------------------------------------------------------------------------------------------"
	print

######################################################################
def get_queued_job_name_list(status_filter, ignore_problem_nodes):

	if(status_filter not in ['ALL', 'QUEUED', 'RUN', 'PEND']): error_exit_with_message("Invalid status_filter (%s)" %(status_filter)) 

	jobs_stat_list = get_queued_jobs_status() #This gives the job_info_list instead of the string.

	job_name_list=[]

	for job_ID in range(len(jobs_stat_list)):

		job_info=copy.deepcopy(jobs_stat_list[job_ID])

		job_state=job_info['STATE']

		job_name=job_info['JOB_NAME']

		if((status_filter=='ALL') or (status_filter=='QUEUED')):
			pass
		elif(status_filter=='RUN'):
			if(job_state!=running_job_state_str()): continue

		elif(status_filter=='PEND'):
			if(job_state!=pending_job_state_str()): continue

		else:
			error_exit_with_message("Invalid status_filter (%s)" %(status_filter))


		######March 09, 2012 IGNORE THE PROBLEM NODES!!###############
		######THESE ARE SPECIFIC TO THE BIOX2-CLUSTER!!###############
		problem_node_list=["node-4-31.local","node-5-13.local", "node-7-1.local", "node-5-36.local"]
	
		if(ignore_problem_nodes==True):

			Is_problem_node=False
			
			#node-6-4.local 1 counts. "node-9-16.local" 1 counts. Will add these to problem node is it clash at least twice.

			#Found  node-4-31 one more time [doesn't write JOB_IB into SLAVE_JOBS/248/JOB_ID.txt]!

			#May 17, 2012. Have problem killing node-5-36.local/7. If this happens again, will add to problem node list.

			for problem_node in problem_node_list:

				if(job_info['EXEC_HOST'].count(problem_node)>0): Is_problem_node=True

			if(Is_problem_node): continue
		##############################################################

		job_name_list.append(job_name)

	return job_name_list


######################################################################
def get_master_job_status(trail_num=1):

	num_trails=0

	master_job_status="BLAH_BLAH_BLAH"

	while(True):

		num_trails+=1

		jobs_stat_list = get_queued_jobs_status() #This gives the job_info_list instead of the string.

		master_tag = abspath( "MASTER" ).replace("/","_")[1:]

		master_job_status="NOT_QUEUED"

		num_master_job_line_found=0

		matching_master_job_line_list=[]

		for job_ID in range(len(jobs_stat_list)):
	
			job_info=copy.deepcopy(jobs_stat_list[job_ID])

			job_state=job_info['STATE']

			job_name=job_info['JOB_NAME']

			if( job_name != master_tag ): continue

			matching_master_job_line_list.append(line_list)

			num_master_job_line_found+=1

			if(job_state==running_job_state_str()): 
				master_job_status="RUNNING"
			elif(job_state==pending_job_state_str()): 
				master_job_status="PENDING"
			else:
				master_job_status="UNKNOWN"	

		if(num_master_job_line_found>1): error_exit_with_message("num_master_job_line_found>1")

		if(master_job_status!="UNKNOWN"): break

		sleep(5)
		if(num_trails>5):
			print "ERROR master job lines:"
			for line_ID in range(len(matching_master_job_line_list)):
				print "matching_master_job_line_list[%d]=" %(line_ID), matching_master_job_line_list[line_ID]

			error_exit_with_message("Unsucessfully determine master_job_status after 5 trails!")

	
	return master_job_status

######################################################################

def get_slave_job_name_list(status_filter):  

	queued_job_name_list=get_queued_job_name_list(status_filter, ignore_problem_nodes=False)

	slave_job_name_list=[]

	generic_slave_tag = abspath( 'SLAVE_JOBS/' ).replace('/','_')[1:] #Assume specific form!

	for n in range( len(queued_job_name_list)):

		job_name=queued_job_name_list[n]

		if(job_name[:len(generic_slave_tag)]==generic_slave_tag):
			slave_job_name_list.append(queued_job_name_list[n])

		if(slave_job_name_list.count(job_name)>1): #Ensure no duplicated job_name
			print_queued_jobs_status()
			print "ERROR: job_name=%s" %(job_name)
			print "ERROR: job_name_list.count(job_name)=%d" %(slave_job_name_list.count(job_name))	
			error_exit_with_message("slave_job_name_list.count(job_name)>1! PERHAPS JOBS FROM PRIOR SUBMISSION WEREN'T PROPERLY REMOVED?")

	return slave_job_name_list


######################################################################
def print_slave_jobs_status():

	jobs_stat_list = get_queued_jobs_status() #This gives the job_info_list instead of the string.

	slave_job_name_list = get_slave_job_name_list('ALL')

	print "---------------------------------------slave_jobs_status:-------------------------------------------"
	print "#%4d: %s" %(0, get_job_stat_colname_line())

	slave_job_count=0

	for job_ID in range(len(jobs_stat_list)):

		job_info=copy.deepcopy(jobs_stat_list[job_ID])

		job_name=job_info['JOB_NAME']

		if(slave_job_name_list.count(job_name)>0):
			slave_job_count+=1
			print "#%4d: %s" %(slave_job_count, job_info_to_job_stat_line(job_info))

	print "-----------------------------------------------------------------------------------------------------"
	print

######################################################################

def print_queued_job_name_list(status_filter):

	job_name_list=get_queued_job_name_list(status_filter, ignore_problem_nodes=False)

	print "---------------------------------------%s_job_name_list------------------------------------------" %(status_filter)

	for n in range(len(job_name_list)):

		print "job_name #%d: %s" %(n+1, job_name_list[n])

	print "-----------------------------------------------------------------------------------------------------"
	print

######################################################################
def print_slave_job_name_list(status_filter):

	job_name_list=get_slave_job_name_list(status_filter)

	print "-------------------------------------%s_slave_job_name_list------------------------------------------" %(status_filter)

	for n in range(len(job_name_list)):

		print "job_name #%d: %s" %(n+1, job_name_list[n])

	print "---------------------------------------------------------------------------------------------------------"
	print


######################################################################
def master_kill_all_slave_jobs_and_exit(exit_message=""):
	
	master_kill_all_slave_jobs_and_exit_scheduler_specific(exit_message)







