#!/usr/bin/env python

from sys import argv,exit
from glob import glob
import string
from os import system,chdir,getcwd
from os.path import basename, dirname, exists, expanduser
from time import sleep
from parse_options import parse_options

# By default remove files.
RM_FILES = not(  parse_options( argv, "no_rm_files", 0 ) )

outfile_cluster = argv[1]
outdir = argv[2]

def wait_for_file( file ):
    for k in range( 10 ):
        if exists( file ):
            break
        else:
            print "waiting for file to show up: ", file
            sleep( 5 )


# Need to be careful to do this WITHIN the directory of interest!
wait_for_file( outfile_cluster )
system( "chmod 777 "+outfile_cluster)

if RM_FILES:
    command = 'rm -rf '+outdir
    print command
    system( command )
