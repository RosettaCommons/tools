#!/usr/bin/env python

from sys import argv,exit
import string,subprocess
from os import system
from os.path import basename,dirname,abspath,exists,expanduser
from cluster_info import *

def Help():
    print
    print argv[0]+' <cluster> [<name of files to sync>]'
    print
    if ( subprocess.call( ['/bin/bash','-i','-c','alias r2c']) == 1 ):
        print "You may want to alias this command to r2c, by putting the following in your .bashrc: "
        print ' alias r2c="rsync_to_cluster.py"'
    exit()

if len(argv)<2:
    Help()

cluster_in = argv[1]
(cluster,remotedir) = cluster_check( cluster_in )

if cluster == 'unknown':
    Help()

args = argv[2:]
dir = ''
extra_args = []

# handle flags like '--exclude' and '--delete' correctly
for m in args:
    if len( m ) > 2 and m.find( '--' ) > -1:
        extra_args.append( m )
    else:
        dir += ' '+m

if len(dir) == 0: dir = '.'

username = basename( expanduser('~') )

# strip off directory name based on local path.
clusterdir = remotedir
clusterdir += strip_home_dirname( abspath('.') )

# make sure directory on cluster is ready for files.
command = 'ssh ' + cluster + ' mkdir -p '+clusterdir
print(command)
system(command)

cluster_prefix = cluster+':'
if len(cluster) == 0: cluster_prefix = ''

# Do it!
command = 'rsync -avzL '+dir+' '+cluster_prefix+clusterdir+' '+string.join(extra_args)
print(command)
system(command)
print
print 'Ran the following command: '
print(command)
