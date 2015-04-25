#!/usr/bin/env python

######################################################################
from scheduler_util import *

######################################################################
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.utility.master_util import *

from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################


#Run this script in the main_folder of the DAG JOB.

JOBDIR = 'SLAVE_JOBS/'
job_prefix = abspath( JOBDIR ).replace('/','_')[1:]

sys.stdout.flush()
sys.stderr.flush()

kill_all_queued_jobs_with_prefix(job_prefix, verbose=True)


