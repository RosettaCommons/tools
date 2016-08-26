#!/usr/bin/env python

from sys import argv,exit
from glob import glob
import string
from os import system
from os.path import basename, dirname, exists, expanduser,abspath
from time import sleep
from parse_options import parse_options

# By default remove files.
RM_FILES = not(  parse_options( argv, "no_rm_files", 0 ) )

indir_prefix = argv[1]
cat_outfile = dirname( indir_prefix) + '/' + basename(indir_prefix).lower() + '_sample.out'
filter_outfile = cat_outfile.replace('.out','.low4000.out')

if exists( filter_outfile ):
    print  argv[0], ' returning early:  Already exists: ', filter_outfile
    exit( 0 )

globstring = indir_prefix+'/*sample.out'
globfiles = glob( globstring )

globstring = indir_prefix+'_S_*/*sample.out'
globfiles2 = glob( globstring )
for file in globfiles2: globfiles.append( file )


#if len( globfiles ) == 0:
#    sleep( 5 )
#    globfiles = glob( globstring )
#globfiles.sort()

if len( globfiles ) == 0:
    print  argv[0], ' returning early:  no files.'
    exit( 0 )

print globfiles


PYDIR = dirname( abspath( argv[0] ) )
assert( exists( PYDIR ) )

command = PYDIR+'/extract_lowscore_decoys_outfile.py '+string.join( globfiles) +' 4000 > '+filter_outfile
print( command )
system( command )
system( 'chmod 777 '+filter_outfile )

if RM_FILES:
    for globfile in globfiles:
        command = 'rm -rf '+globfile
        print( command )
        system( command )


