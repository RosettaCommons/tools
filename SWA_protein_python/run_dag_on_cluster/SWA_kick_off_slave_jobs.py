#!/usr/bin/env python

from sys import argv,exit,stdout
from SWA_dagman_LSF_continuous import *
from SWA_util import *

try:
    N_JOBS = int( argv[ 1 ] )
except:
    print 'must supply argument:\n   %s <number slave jobs>' % argv[0]
    exit()

job_cluster_number = kick_off_slave_jobs( N_JOBS )
