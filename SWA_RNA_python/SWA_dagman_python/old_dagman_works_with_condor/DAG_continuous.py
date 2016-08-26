#!/usr/bin/env python

from sys import argv,exit,stdout
import string
from os.path import basename,exists
from random import randrange
from time import sleep
######################################################################

from DAG_continuous_util import *

######################################################################
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################

# In this implementation, kick off a certain number
#  of long-running scripts -- they will carry out communication
#  this master script based on files that show up in
#  job-specific directories.

START_argv=copy.deepcopy(argv)

print "ENTER: %s" %(list_to_string(START_argv))

do_delete_log_files=parse_options( argv, "do_delete_log_files", "True" )

N_JOBS=parse_options( argv, "j", 100)

print "do_delete_log_files=%s " %(do_delete_log_files)

if(len(argv)!=2): error_exit_with_message("len(argv)!=2, leftover_argv=%s" %(list_to_string(argv) ) )

dagman_file = argv[1]

#########################################
# Parse dagman file
#########################################
try:
	lines = open( dagman_file ).readlines()
except:
	kill_all_slave_jobs_and_exit("Cannot open dagman_file: %s" %(dagman_file))

jobs = [] #list of jobs
condor_submit_file = {} #map from job to the condor filename
post_script = {} #map from job to post_script command
pre_script = {} #map from job to pre_script command
parents = {} #map from child to a list of parent
for line in lines:
    print line[:-1]
    stdout.flush()
    if len( line ) > 4 and line[:3] == "JOB":
        cols = string.split( line )
        job = cols[1]
        jobs.append( job )
        condor_submit_file[ job ] = cols[2]
    elif len( line ) > 6 and line[:6] == "SCRIPT":
        cols = string.split( line )
        job = cols[2]
        if cols[1] == "PRE":
            pre_script[ job ] = string.join( cols[3:] )
        else:
            assert( cols[1] == "POST" )
            post_script[ job ] = string.join( cols[3:] )
    elif len( line ) > 7 and line[:6] == "PARENT":
        cols = string.split( line )
        assert( cols.count( "CHILD" ) )
        child_index = cols.index( "CHILD" )
        parent_jobs =  cols[1:child_index]
        print "child_index " , child_index , " parent_jobs =" , parent_jobs
        child_jobs =  cols[child_index:]
        for child_job in child_jobs:
            if child_job not in parents.keys():
                parents[ child_job ] = []
            for parent_job in parent_jobs: parents[ child_job ].append( parent_job ) 
            
#if a job has not previously been indicated to be a child, it cannot have a parent job
#however this should be ok since a job always have itself as one of its child            
#map[key]=value

# Kick off jobs!
job_cluster_number = kick_off_slave_jobs( N_JOBS )
#sleep( 10 ) #OK give time for each node to check if it is broken

##check that the job_names in jobs are unique
prev_jobs=[]
for n in range(len(jobs)):
	if(jobs[n] in prev_jobs): kill_all_slave_jobs_and_exit("jobs[n] (%s) already exist in prev_jobs!" %(jobs[n]) )
	prev_jobs.append(jobs[n])

prev_jobs[:]=[] #clear

####CLUSTERING JOBS have higher priority, so MOVE ALL CLUSTERING JOBS TO THE TOP OF THE LIST###

old_jobs=jobs
jobs= []

for job in old_jobs: #cluster jobs

	if(job in jobs): continue #already added!

	if( (condor_submit_file[job].count('cluster')>0) or (condor_submit_file[job].count('CLUSTER')>0) ): jobs.append( job)

for job in old_jobs: #A VIRT_RIBOSE_SAMPLER job

	if(job in jobs): continue #already added!

	if(job.count('VIRT_RIBOSE_SAMPLER')>0): jobs.append( job)


for job in old_jobs: #A filterer job.

	if(job in jobs): continue #already added!

	if(job.count('FILTERER')>0): jobs.append( job)


for job in old_jobs: #other jobs

	if(job in jobs): continue #already added!

	jobs.append( job)

######debug check########################################################
if( len(jobs)!=len(old_jobs) ):
	kill_all_slave_jobs_and_exit("len(jobs)!=len(old_jobs)")
		
for job in jobs:
	if(old_jobs.count(job)!=1):
		kill_all_slave_jobs_and_exit("old_jobs.count(%s)!=1" %(job))

old_jobs[:]=[] #Clear
 
print "Print jobs for debugging"		
for job in jobs:
	print condor_submit_file[job]
print "------------------------"

###############################################################################################
print_title_text("Start job submissions")

done = {} #map of job to whether it is done: 0 not-done/false, 1 done/true
queued = {} #map of job to wether it has been queued yet: 0 false, 1 true 
for job in jobs:
    done[ job ] = False
    queued[ job ] = False

check_slave_node_disappear_wait_count=0

job_info_list_map = {} #A map from job to the out_log_files of the job
all_done = False
while not all_done:  

	###################################################
	# Find jobs that are ready to go but not done.
	###################################################
	all_done = True 
	for job in jobs: 
#		print condor_submit_file[job]
#		sys.stdout.flush()
#		sys.stderr.flush()

		queued_a_job = False
		if not done[ job ]:
			all_done = False 

			if not queued[ job ]: #Consider queuing the job

				ready_to_queue = True 
 				if job in parents.keys():  #if job is a child
					for parent_job in parents[job]:  #iterate over the job's parent
						if not done[ parent_job ]:   #Check if parent job is done
							ready_to_queue = False
							break

				if ready_to_queue:
					if job in pre_script.keys():
						pre_command = pre_script[ job ]
						pre_command_log_file = condor_submit_file[job] +'.pre_script.log'
						command =  pre_command + ' > '+pre_command_log_file  #delete any old pre_command_log_file in the process
						submit_subprocess( command, True )
			
					post_script_command=""
					if job in post_script.keys():
						post_script_command = post_script[ job ] 
		
					job_info_list_map[ job ] = condor_submit( condor_submit_file[ job ], post_script_command )
					queued[ job ] = True
 					break # some other jobs may have finished while we were queuing. TURN THIS ON, ...AND start for loop from 1st job again....
 				

	sleep(1)

	if(check_slave_node_disappear_wait_count==300): check_slave_node_disappear_wait_count=0
	check_slave_node_disappear_wait_count+=1

  ###################################################
  # Check for any jobs that are done
  ###################################################
	error_check_all_slaves(verbose=False)

	for job in jobs:
		if not done[ job ] and queued[ job ]:

			job_info_list=job_info_list_map[ job ]

			if(check_slave_node_disappear_wait_count==300): #as to not call bjobs command too frequently
				check_if_slave_node_disappear(job_info_list)

#			if(len(job_info_list)==0):
#				print "job= " ,job, " job_info_list= ", job_info_list
#				kill_all_slave_jobs_and_exit("len(job_info_list)==0 for job %s " %(job) )

			job_is_done = check_if_job_in_done( job_info_list )
	
			stdout.flush()

			if(job_is_done):
				done[ job ] = True
				queued[ job ] = False

				if(do_delete_log_files): delete_log_files(job_info_list, keep_some_log_file=True)


send_finish_signal_to_slave_nodes(N_JOBS)


print "Successfuly RAN DAG_continuous.py!!"
