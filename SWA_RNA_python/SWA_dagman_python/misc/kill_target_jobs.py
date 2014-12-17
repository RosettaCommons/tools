#!/usr/bin/python

from sys import argv
from os import popen

assert( len(argv) > 1 )
target = argv[1]

#death_row = {}
death_row = []
qstat = popen( 'qstat' ).readlines()

for line in qstat:
    if '.' not in line: continue
    job_id = line.split('.')[0]
    full_job_info = popen( 'qstat -f '+job_id ).read()
    if target not in full_job_info: continue
    death_row.append( job_id )

for job_id in death_row:
    executed = popen( 'qdel '+job_id ).read()
    print 'executed job_id: '+job_id


