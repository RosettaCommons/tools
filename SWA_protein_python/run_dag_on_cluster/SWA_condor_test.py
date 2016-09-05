#!/usr/bin/env python

from sys import argv,exit,stdout
from os import system
from SWA_dagman_LSF_continuous import *

condor_submit_file_ = argv[1]

print_title_text('Run: '+ condor_submit_file_)
sys.stdout.flush()
sys.stderr.flush()

lines = safe_open( condor_submit_file_ ).readlines()
log = ""
output = ""
err = "/dev/null"
exe = ""
args = ""
universe = "vanilla"
queue_num = 0
for line in lines:
	if len( line ) > 2:
		cols = string.split( line )
		if cols[0] ==  "executable":
			assert( cols[1] == "=" )
			exe = cols[2]
		elif cols[0] == "arguments":
			assert( cols[1] == "=" )
			args = string.join(cols[2:])
		elif cols[0] == "log":
			assert( cols[1] == "=" )
			log = cols[2]
		elif cols[0] == "output":
			assert( cols[1] == "=" )
			output = cols[2]
		elif cols[0] == "error":
			assert( cols[1] == "=" )
			err = cols[2]
		elif cols[0] == "universe":
			assert( cols[1] == "=" )
			universe = cols[2]
		elif (cols[0]).lower() == "queue":
			if len( cols ) > 1: queue_num = int( cols[1] )

q_string = '%d' % 0
args_new = args.replace('$(Process)',q_string).replace('$(PROCESS)',q_string)
err_new = err.replace('$(Process)',q_string).replace('$(PROCESS)',q_string)
output_new = output.replace('$(Process)',q_string).replace('$(PROCESS)',q_string)

command = "%s %s"  % (exe, args_new )
print command
system( command )

