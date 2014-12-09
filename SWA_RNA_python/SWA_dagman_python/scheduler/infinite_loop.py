#!/usr/bin/env python

######################################################################
from scheduler_util import *

######################################################################
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.utility.master_util import *

from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################



count=0

while(True):
	count+=1

	print "count=%d" %(count) 
	sys.stdout.flush()
	sys.stderr.flush()
	sleep(100)

