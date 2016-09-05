#!/usr/bin/env python

from sys import argv,exit,stdout
import string
from os.path import basename,exists
from random import randrange
from time import sleep
from parse_options import parse_options

from SWA_dagman_LSF_continuous import *
from SWA_util import *

# In this implementation, kick off a certain number
#  of long-running scripts -- they will carry out communication
#  this master script based on files that show up in
#  job-specific directories.
N_JOBS = parse_options( argv, "j", 100 );
finalize = not parse_options(argv, "no_finalize", 0 )

# This option determines the point at which we queue more slave jobs (up to N_JOBS).
num_start_slaves = parse_options( argv, "num_start_slaves", 1 )
all_slaves_at_once = parse_options( argv, "all_slaves_at_once", 0 )

nhours_slave = parse_options( argv, "nhours", 48 )

#########################################
# Parse dagman file
#########################################
dagman_file = argv[1]
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
    #print line[:-1]
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
        child_jobs =  cols[child_index:]
        for child_job in child_jobs:
            if child_job not in parents.keys():
                parents[ child_job ] = []
            for parent_job in parent_jobs: parents[ child_job ].append( parent_job )

#if a job has not previously been indicated to be a child, it cannot have a parent job
#however this should be ok since a job always have itself as one of its child
#map[key]=value

# Kick off jobs!
if ( all_slaves_at_once ):
    job_cluster_number = kick_off_slave_jobs( N_JOBS, nhours_slave )
else:
    job_cluster_number = kick_off_slave_jobs( min( N_JOBS, num_start_slaves  ), nhours_slave  )


sleep( 5 ) # Wait for some of the jobs to register!!!

####CLUSTERING JOBS have higher priority, so MOVE ALL CLUSTERING JOBS TO THE TOP OF THE LIST###
old_jobs=jobs
jobs= []

for job in old_jobs: #cluster jobs
	if(condor_submit_file[job].count('cluster')>0):
		jobs.append( job)

for job in old_jobs: #other jobs
	if(condor_submit_file[job].count('cluster')==0):
		jobs.append( job)

######debug check########################################################
if( len(jobs)!=len(old_jobs) ):
	kill_all_slave_jobs_and_exit("len(jobs)!=len(old_jobs)")

for job in jobs:
	if(old_jobs.count(job)!=1):
		kill_all_slave_jobs_and_exit("old_jobs.count(%s)!=1" %(job))

#print "Print jobs for debugging"
#for job in jobs:
#	print condor_submit_file[job]
###############################################################################################
print_title_text("Start job submissions")

done = {} #map of job to whether it is done: 0 not-done/false, 1 done/true
queued = {} #map of job to wether it has been queued yet: 0 false, 1 true
for job in jobs:
    done[ job ] = False
    queued[ job ] = False

output_files = {} #A map from job to the output_files of the job
all_done = False
num_active_jobs = 0
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
						if not exists( pre_command_log_file ):
							submit_subprocess( command, True )

					output_files[ job ] = condor_submit( condor_submit_file[ job ], N_JOBS, nhours_slave )
					queued[ job ] = True
 					break # some other jobs may have finished while we were queuing. TURN THIS ON, ...AND start for loop from 1st job again....


	sleep(1)

    ###################################################
    # Check for any jobs that are done
    ###################################################
        num_active_jobs = 0
	for job in jobs:
		if not done[ job ] and queued[ job ]:
			still_running = check_output_files( output_files[ job ] )

			if len(output_files[job]) > 0 : print "Jobs still running: ",output_files[job][-1], still_running
			else: print 'No outfiles: ', job

			stdout.flush()

			if not still_running:
				if job in post_script.keys(): # and len( output_files[job] ) > 0:   #The last clause says don't run postscript if the job was already done!
					command = post_script[ job ]
					print( command )
					submit_subprocess( command, True )
				done[ job ] = True
				queued[ job ] = False
                        else:
                            num_active_jobs += 1

        if (num_active_jobs > 1): print "----------------------------------------------"

if finalize: finish_jobs()

