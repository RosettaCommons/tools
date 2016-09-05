#!/usr/bin/env python

from sys import argv, exit, stdout
import string 
from os.path import basename, exists
from random import randrange
from time import sleep
###############################################################################
from DAG_continuous_qsub_util import *
###############################################################################
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
###############################################################################

# In this implementation, kick off a certain number of long-running 
# scripts -- they will communicate with this master script based on
# files that show up in job-specific directories

START_argv=copy.deepcopy(argv)

print "ENTER: %s" % list_to_string(START_argv)

do_delete_log_files=parse_options( argv, "do_delete_log_files", "True" )
N_JOBS=parse_options( argv, "j", 100 )

print "do_delete_log_files=%s " % do_delete_log_files

if(len(argv)!=2): error_exit_with_message("len(argv)!=2, leftover_argv=%s" % list_to_string(argv))

dagman_file = argv[1]

###############################################################################
# Parse dagman file
###############################################################################
try:
	lines = open( dagman_file ).readlines()
except:
	kill_all_slave_jobs_and_exit("Cannot open dagman_file: %s" % dagman_file)

jobs = []
condor_submit_file = {}
post_script = {}
pre_script = {}
parents = {}
for line in lines:
	print line[:-1]
	stdout.flush()
	cols = line.split(' ')
	if cols[0] == "JOB":
		job = cols[1]
		jobs.append( job )
		condor_submit_file[ job ] = cols[2]
	elif cols[0] == "SCRIPT":
		job = cols[2]
		if cols[1] == "PRE":
			pre_script[ job ] = string.join( cols[3:] )
		else:
			assert( cols[1] == "POST" )
			post_script[ job ] = string.join( cols[3:] )
	elif cols[0] == "PARENT":
		assert( cols.count( "CHILD" ) )
		child_index = cols.index( "CHILD" )
		parent_jobs = cols[1:child_index]
		print "child_index ", child_index, " parent_jobs = ", parent_jobs
		child_jobs = cols[child_index:]
		for child_job in child_jobs:
			if child_job not in parents.keys():
				parents[ child_job ] = []
			for parent_job in parent_jobs:
				parents[ child_job ].append( parent_job )

###############################################################################
# Kick off jobs!
###############################################################################
job_cluster_number = kick_off_slave_jobs( N_JOBS )

### Check that the job_names in jobs are unique
prev_jobs = []
for n in range(len(jobs)):
	if(jobs[n] in prev_jobs): kill_all_slave_jobs_and_exit("jobs[n] (%s) already exist in prev_jobs!" %(jobs[n]) )
	prev_jobs.append( jobs[n] )

### Clear prev_jobs
prev_jobs[:] = []

### Move all clustering jobs to the top of the list, they have highest priority
old_jobs = jobs
jobs = []

# CLUSTERER jobs
for job in old_jobs:
	if( job in jobs ): continue
	if( condor_submit_file[job].count('cluster') or condor_submit_file[job].count('CLUSTER') ): jobs.append( job)

# VIRT_RIBOSE_SAMPLER jobs
for job in old_jobs:
	if( job in jobs ): continue
	if( job.count('VIRT_RIBOSE_SAMPLER') ): jobs.append( job )

# FILTERER jobs
for job in old_jobs:
	if( job in jobs ): continue
	if( job.count('FILTERER') ): jobs.append( job )

# OTHER jobs
for job in old_jobs:
	if( job in jobs ): continue
	jobs.append( job )


######debug check########################################################
if( len(jobs)!=len(old_jobs) ):
	kill_all_slave_jobs_and_exit("len(jobs)!=len(old_jobs)")
		
for job in jobs:
	if(old_jobs.count(job)!=1):
		kill_all_slave_jobs_and_exit("old_jobs.count(%s)!=1" %(job))

# Clear old_jobs
old_jobs[:]=[]
 
print "Print jobs for debugging"		
for job in jobs:
	print condor_submit_file[job]
print "------------------------"


###############################################################################
# Start Job Submissions
###############################################################################
print_title_text("Start job submissions")

done = {}
queued = {}
for job in jobs:
	done[ job ] = False
	queued[ job ] = False

check_slave_node_disapear_wait_count=0
job_info_list_map = {}
all_done = False

while not all_done:

	###########################################################################
	# Find jobs that are ready to go but not done
	###########################################################################
	all_done = True
	for job in jobs:
		
		queued_a_job = False
		
		if not done[ job ]:
		
			all_done = False
		
			if not queued[ job ]:
			
				ready_to_queue = True
		
				if job in parent.keys():
					for parent_job in parents[ job ]:
						if not done[ parent_job ]:
							ready_to_queue = False
							break

				if ready_to_queue:
					if job in pre_script.keys():
						pre_command = pre_script[ job ]
						pre_command_log_file = condor_submit_file[ job ] + '.pre_script.log'
						command = pre_command + ' > ' + pre_command_log_file
						submit_subprocess( command, True )

					post_script_command=''
					if job in post_script.keys():
						post_script_command = post_script[ job ]

					job_info_list_map[ job ] = condor_submit( condor_submit_file[ job ], post_script_command )
					queued[ job ] = True
					break


	sleep(1)

	if( check_slave_node_disapear_wait_count==300 ):	check_slave_node_disapear_wait_count=0
	check_slave_node_disappear_wait_count+=1

	###########################################################################
	# Find jobs that are ready to go but not done
	###########################################################################
	error_check_all_slaves(verbose=False)

	for job in jobs:
		if not done[ job ] and queued[ job ]:

			job_info_list = job_info_list_map[ job ]

			if( check_slave_node_disappear_wair_count == 300 ):
				check_if_slave_node_disappear( job_info_list )

			job_is_done = check_if_job_in_done( job_info_list )
			stdout.flush()

			if( job_is_done ):
				done[ job ] = True
				queued[ job ] = False

				if( do_delete_log_files ): delete_log_files( job_info_list, keep_some_log_file=True )


send_finish_signal_to_slave_nodes(N_JOBS)
print "Successfully RAN DAG_continuous_qsub.py!!"				
































