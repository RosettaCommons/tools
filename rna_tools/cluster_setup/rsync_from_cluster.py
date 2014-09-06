#!/usr/bin/python

from sys import argv,exit
import string,subprocess
from os import system
from os.path import basename,dirname,abspath,exists,expanduser
from cluster_info import *

def Help():
    print
    print argv[0]+' <cluster> <any extra rsync flags>'
    print
    if ( subprocess.call( ['/bin/bash','-i','-c','alias rfc']) == 1 ):
        print "You may want to alias this command to rfc, by putting the following in your .bashrc: "
        print ' alias rfc="rsync_from_cluster.py"'
    exit()

if len(argv)<2:
    Help()

cluster_in = argv[1]
(cluster,remotedir) = cluster_check( cluster_in )

if cluster == 'unknown':
    Help()

# handle flags like '--exclude' and '--delete' correctly
args = argv[2:]
filename = ''
extra_args = []
for m in args:
    if len( m ) > 2 and m.find( '--' ) > -1:
        extra_args.append( m )
    else:
        assert( len(filename) == 0 )
        filename = m

dir = '.'

clusterdir = abspath(dir).replace('/Users/%s/' % user_name,'')
clusterdir = clusterdir.replace('/scratch/users/%s/' % user_name,'')
clusterdir = clusterdir.replace('/work/%s/' % user_name,'')
clusterdir = clusterdir.replace('/home1/%s/%s/' % ( xsede_dir_number, xsede_user_name ),'')
clusterdir = clusterdir.replace('/work/%s/%s/' % (xsede_dir_number, xsede_user_name ),'')

clusterdir = remotedir+clusterdir

cluster_prefix = cluster+':'
if len(cluster) == 0: cluster_prefix = ''

command = 'rsync -avzL '+cluster_prefix+clusterdir + '/' + filename + ' . '+string.join(extra_args)+' --exclude="condor*log"'
print(command)
system(command)

