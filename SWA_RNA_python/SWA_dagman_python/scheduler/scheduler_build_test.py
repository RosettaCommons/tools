#!/usr/bin/env python

######################################################################
from SWA_dagman_python.utility.SWA_util import *

from SWA_dagman_python.dagman.DAG_continuous_util import *
######################################################################



#Use this script to check whether the scheduler scripts are working properly!

#Note our scheduler script keep tracks and manage jobs mainly through the job_name variable which is set to equal to:
#SLAVE_DIR = 'SLAVE_JOBS/%d' %(job_num) where job_num = 0,1,2,3.., N_JOBS-1
#slave_job_name = abspath( SLAVE_DIR ).replace('/','_')[1:]
#So if run this script at /home/john/ then slave_job_name= home_john_SLAVE_JOBS_0, home_john_SLAVE_JOBS_1, .... home_john_SLAVE_JOBS_{N_JOBS-1}

#If you ran this script and it failed, make sure to wait at least 2 minutes
#before rerunning the script (this give enough time for all prior submitted jobs to run to completion)
#OR you can just run this script in another directory.
######################################################################

#~/SWA_RNA_python/SWA_dagman_python/scheduler/scheduler_build_test.py

print "Enter scheduler_build_test.py"
print
######1. OK try to submit 5 slave jobs########## 
print "CHECKING JOB SUBMISSION COMMANDS:"


N_JOBS=5

if(exists('SLAVE_JOBS')): submit_subprocess("rm -r SLAVE_JOBS/")

kick_off_slave_jobs( N_JOBS, SLAVE_EXE = get_PYEXE('scheduler/simple_test_script.py'), wall_time=1)

sleep(2)
######2. Check that the jobs exist!########## 
print "CHECKING JOB STATUS COMMANDS:"
print 

print_slave_jobs_status()

print_slave_job_name_list("QUEUED")

print_slave_job_name_list("RUN")

print_slave_job_name_list("PEND")

queued_slave_job_list= get_slave_job_name_list("QUEUED")

running_slave_job_list=get_slave_job_name_list("RUN")

pending_slave_job_list=get_slave_job_name_list("PEND")

num_slave_job_queued=len(queued_slave_job_list)

num_slave_job_running=len(running_slave_job_list)

num_slave_job_pending=len(pending_slave_job_list)

sleep(5) #Give sometime for the simple_test_script.py to run!
######3. Kill submitted jobs!########## 
print "CHECKING JOB KILLING COMMANDS:"
print 

JOBDIR = 'SLAVE_JOBS/'
job_prefix = abspath( JOBDIR ).replace('/','_')[1:]

num_slave_job_killed=kill_all_queued_jobs_with_prefix(job_prefix, verbose=True)

sleep(20)
######4. Print results
print "-----------------------------------SCHEDULER_BUILD_TEST_RESULTS:-----------------------------------------------"
print

if(num_slave_job_queued!=N_JOBS):
	print "WARNING: Script queued ONLY %d of %d JOBS!" %(num_slave_job_queued, N_JOBS)
else:
	print "SUCCESS: Script queued all %d of %d JOBS!" %(num_slave_job_queued, N_JOBS)

print

if(num_slave_job_running!=N_JOBS):
	print "WARNING: ONLY %d of %d JOBS were running!" %(num_slave_job_running, N_JOBS)
else:
	print "SUCCESS: All %d of %d JOBS were running!" %(num_slave_job_running, N_JOBS)

print

if(num_slave_job_pending!=0):
	print "WARNING: %d of %d JOBS were pending!" %(num_slave_job_pending, N_JOBS)
else:
	print "SUCCESS: %d of %d JOBS were pending!" %(num_slave_job_pending, N_JOBS)

print

if(num_slave_job_killed!=N_JOBS):
	print "WARNING: Script killed ONLY %d of %d JOBS!" %(num_slave_job_killed, N_JOBS)
else:
	print "SUCCESS: Script killed %d of %d JOBS!" %(num_slave_job_killed, N_JOBS)

print

after_kill_slave_job_name_list= get_slave_job_name_list("ALL")

if(len(after_kill_slave_job_name_list)>0): 
	print "WARNING: After killing all slave_jobs, %s SLAVE_JOBS still remains in the queue!"  %( len(after_kill_slave_job_name_list) )
	print
	print_slave_jobs_status()

print "-----------------------------------------------------------------------------------------------------"
print "EXIT scheduler_build_test.py"
