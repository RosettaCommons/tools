#!/usr/bin/env python

from sys import argv
from os import system

def Help():
    print
    print argv[0], ' <Puzzle Number> '
    print
    print '  pdb file names should be specified in order in submit_order.list'
    print
    exit( 0 )
if len( argv ) == 1:
    Help()


lines = open( 'submit_order.list' ).readlines()
pdbs = map( lambda x:x.replace( '\n', ''), lines )

problem_number = int( argv[1] )

count = 0
for pdb in pdbs:
    count = count+1

    new_pdb = 'DASLAB_Problem%d_Rank%d.pdb' % (problem_number, count )
    command = 'reorder_to_standard_pdb.py %s > %s' % ( pdb, new_pdb )
    print command
    system( command )

command = 'tar cvfz DASLAB_Problem%d.tgz DASLAB_Problem%d_Rank*.pdb' % (problem_number, problem_number)
print command
system( command )
