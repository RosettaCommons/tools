#!/usr/bin/python

from sys import argv, exit
from os import popen
import subprocess

MAX_TRIES = 3

assert( len(argv) > 1 )
target = argv[1]
safe_mode = False
if len(argv) > 2:
    safe_mode = ('safe_mode' in argv[2])
death_row = []

count = 0
while True:
    qstat_full_info, err = subprocess.Popen( ['qstat','-f'], stdout=subprocess.PIPE,stderr=subprocess.PIPE ).communicate()
    if (count > MAX_TRIES) or (not err) or (not len(err)): break
    count += 1

qstat_full_info = qstat_full_info.split('\n\n')
for full_job_info in qstat_full_info:
    if target not in full_job_info: continue
    full_info = full_job_info.split('\n')
    for i, info_line in enumerate( full_info ):
        if 'Job Id' in info_line:
            job_id = info_line.split(':')[1].replace(' ','')
            continue
        if 'Job_Name' not in info_line: continue
        job_name = info_line.split('=')[1]
        for j in xrange(1, 5):
            if '=' in full_info[i+j]: break
            job_name+= full_info[i+j]
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

found_slaves = len(death_row) - int(found_master)

print 'found_slaves: ', found_slaves
print 'found_master: ',found_master

if not safe_mode:
    while True:#not found_master:
        kill = raw_input('Would you like to kill jobs [y/N]: ').lower()
        if 'n' in kill:
            break
        if 'y' in kill:
            for ( job_id, job_name ) in death_row:
                count = 0
                while True:
                    executed, err = subprocess.Popen( ['qdel',job_id], stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()
                    if (count > MAX_TRIES) or (not err) or (not len(err)): break
                    count += 1
            break
