#!/usr/bin/env python

######################################################################
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.scheduler.scheduler_util import *

from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################


START_argv=copy.deepcopy(argv)

print "---------------------------------------------------------------------------------------------------------------------"
print "---------------------------------------------------------------------------------------------------------------------"
print "Enter %s" %(list_to_string(START_argv)) 
print

print "SCHEDULER_TYPE=%s" %(SCHEDULER_TYPE)

#/Users/sripakpa/SWA_RNA_python/SWA_dagman_python/dagman/submit_DAG_job.py -master_wall_time 48 -master_memory_reserve 2048 -num_slave_nodes 50

master_wall_time  = parse_options( argv, "master_wall_time", 48 ) #In hours

master_memory_reserve = parse_options( argv, "master_memory_reserve", 0 ) 

num_slave_nodes = parse_options( argv, "num_slave_nodes", 0 ) 

dagman_file = parse_options( argv, "dagman_file", "" ) 

verbose = parse_options( argv, "verbose", "True" ) 

submit_DAG_job_scheduler_specific(master_wall_time, master_memory_reserve, num_slave_nodes, dagman_file, verbose)

print "------------------------------------------------------------------------------------------"
print "SUCCESSFULLY RAN schedule_util.py!"
print "------------------------------------------------------------------------------------------"
