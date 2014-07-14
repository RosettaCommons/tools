#!/usr/bin/python

from sys import argv
from os import system

lines = open( 'submit_order.list' ).readlines()
pdbs = map( lambda x:x[:-1], lines )

count = 0
for pdb in pdbs:
    count = count+1

    command = 'rna_rosetta_to_pdb.py  -pdb %s ' % pdb
    print command
    system( command )

    new_pdb = pdb.replace( '.pdb', '_standard.pdb' )
    command = 'mv %s DASLAB_%04d.pdb' % (new_pdb, count )
    print command
    system( command )
