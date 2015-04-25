#!/usr/bin/env python

######################################################################
from scheduler_util import *

######################################################################
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.utility.master_util import *

from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################


START_argv=copy.deepcopy(argv)

print "---------------------------------------------------------------------------------------------------------------------"
print "---------------------------------------------------------------------------------------------------------------------"
print "Enter %s" %(list_to_string(START_argv)) 
print

#/Users/sripakpa/SWA_RNA_python/SWA_dagman_python/scheduler/scheduler_queue_job.py -wall_time 144 -outfile master_log.out -errfile master_log.err -job_name Users_sripakpa_minirosetta_test_March_31_SWA_minimize_QSUB_test_trail_1_main_folder_SWA_denovo_CUUG_tetraloop_MASTER /Users/sripakpa/SWA_RNA_python/SWA_dagman_python/dagman/DAG_continuous.py -j 50 rna_build.dag

wall_time  = parse_options( argv, "wall_time", 48 ) #In hours

outfile = parse_options( argv, "outfile", "default_log.out" )

errfile = parse_options( argv, "errfile", "default_log.err" ) 

job_name = parse_options( argv, "job_name", "default_job_name" ) 

job_dir_name = parse_options( argv, "job_dir_name", "" )  #Relative path to the directory the python script is executable. If job_dir_name=="", then job is executed in the current directory.

memory_limit = parse_options( argv, "memory_limit", 4096 ) 

memory_reserve = parse_options( argv, "memory_reserve", 0 ) 

queue_name = parse_options( argv, "queue_name", "DEFAULT" ) 

verbose = parse_options( argv, "verbose", "True" ) 

job_script = list_to_string(argv[1:])

queue_job(job_name, outfile, errfile, job_script, job_dir_name, wall_time, memory_limit, memory_reserve, queue_name, verbose)


print "------------------------------------------------------------------------------------------"
print "SUCCESSFULLY RAN scheduler_queue_job.py!"
print "------------------------------------------------------------------------------------------"
