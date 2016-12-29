step#!/usr/bin/env python

import string
from sys import argv
from os import system,popen
from os.path import exists,basename

outdir = argv[1]
dir_prev = argv[2]
condor_submit_file = argv[3]

if len( argv ) > 4:
    sub_job_tag = argv[ 4 ]
else:
    sub_job_tag = 'START_FROM_'+dir_prev.upper()

# Need directories to be setup. Could do this with a preprocessing script? That would also allow
# for convenient setup of Queue number in the condor script.
MAX_JOBS = 4000
lines = popen( 'grep SCORE %s_sample.cluster.out' % ( dir_prev.lower() ) ).readlines()

tags = map( lambda x: string.split( x )[-1], lines )

for q in range( MAX_JOBS ):
    tag = 'S_%d' % q

    if tag in tags:
        newdir = outdir+'/%s_%s' % ( sub_job_tag, tag )
        if not exists( newdir ):
            system( 'mkdir -p '+newdir )
            system( 'chmod 777 '  + newdir )
    else:
        break

N_JOBS = q

# Go through condor submission file and update number of jobs...
# May need to be careful about any special cases.
lines = open( condor_submit_file ).readlines()
fid = open( condor_submit_file, 'w' )

for line in lines:
    if len( line ) > 5 and line[:5] == 'Queue':
        fid.write( 'Queue %d\n' % N_JOBS )
    else:
        fid.write( line )

fid.close()

