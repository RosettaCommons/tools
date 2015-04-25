#!/usr/bin/env python

######################################################################
from scheduler_util import *

######################################################################
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.utility.master_util import *

from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################

HOMEDIR = expanduser('~')

os.chdir( HOMEDIR )

print_queued_jobs_status()

