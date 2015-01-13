#!/usr/bin/python

from sys import argv, exit
from os import popen

assert( len(argv) > 1 )
target = argv[1]
death_row = []
qstat = popen( 'qstat' ).readlines()

for line in qstat:
    if '.' not in line: continue
    job_id = line.split(' ')[0]
    full_job_info = popen( 'qstat -f '+job_id ).read()
    if target not in full_job_info: continue
    full_job_info = full_job_info.split('\n')
    for i, info_line in enumerate(full_job_info):
        if 'Job_Name' not in info_line: continue
        job_name = info_line.split('=')[1]
        if '=' not in full_job_info[i+1]:
            job_name += full_job_info[i+1]
        if '=' not in full_job_info[i+2]:
            job_name += full_job_info[i+2]
    job_name = job_name.replace('\n','').replace('\t','').replace(' ','')
    job_name = job_name[job_name.index(target):]
    death_row.append( (job_id, job_name) )

found_target_jobs = len( death_row )
assert( found_target_jobs )

found_master = False
for ( job_id, job_name ) in death_row:
    if 'MASTER' in job_name:
        found_master = True
    print '  job_id: ',job_id
    print 'job_name: ',job_name

print 'found_master: ',found_master

while True:#not found_master:
    kill = raw_input('Would you like to kill jobs [y/N]: ').lower()
    if 'n' in kill:
        break
    if 'y' in kill:
        for ( job_id, job_name ) in death_row:
            executed = popen( 'qdel '+job_id ).read()
        break
