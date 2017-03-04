#!/usr/bin/python
import subprocess
import glob
from os.path import exists

queue = {}
condor_files = glob.glob( "CONDOR/SAMPLER/*.condor" )
for condor_file in condor_files:
    cmd = ['grep', 'Queue', condor_file]
    out,err = subprocess.Popen( cmd,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE).communicate()
    key = condor_file.split('/')[-1].split('.')[0]
    queue[ key ] = int(out.replace('\n','').split(' ')[-1])

print '\nQUEUE:'
for key in sorted(queue.keys()):
    start_from_dir = key.replace( '_START', '/START' )
    done_signal = ( '%s/deleting_files_signal_%s.txt' %
                    ( start_from_dir.split('/')[0],
                      start_from_dir.split('/')[-1].lower() ) )
    if not exists( start_from_dir ) and exists( done_signal ):
        queue.pop(key)
        continue
    print key+': ', queue[ key ]

print '\nJOBS COMPLETED:'
for key in sorted(queue.keys()):
    outfiles = glob.glob( key.replace( '_START', '/START' ) + '/*/*out' )
    noutfiles = len( outfiles )
    print key+':', noutfiles, 'of', queue[ key ]
    queue[ key ] -= noutfiles

print '\nJOBS REMAINING:'
total = 0
for key in sorted(queue.keys()):
    print key+':', queue[ key ]
    total += queue[ key ]

print '\nTOTAL JOBS LEFT:', total

#regions_left = glob( 'REGION*/START*' )
#for region in regions_left:
#    print region
