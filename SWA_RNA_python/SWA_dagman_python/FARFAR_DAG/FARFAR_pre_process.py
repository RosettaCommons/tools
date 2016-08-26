#!/usr/bin/env python

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
import os
from os.path import basename, dirname, exists, expanduser
from time import sleep
from sets import Set
import copy
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
from SWA_dagman_python.utility.DAGMAN_util import update_CONDOR_file_with_actual_job_queue_num
######################################################################


python_command=list_to_string(argv)
print_title_text("Enter: " + python_command)


num_node_per_cutpoint = int(argv[1]) 
condor_submit_file = argv[2]

update_CONDOR_file_with_actual_job_queue_num(condor_submit_file, num_node_per_cutpoint)


print_title_text("Exit: " + python_command)
