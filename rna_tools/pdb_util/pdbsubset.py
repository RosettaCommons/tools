#!/usr/bin/python

from os import popen,system
from os.path import basename,exists
from sys import argv,stderr
import string
from read_pdb import read_pdb
from parse_options import get_resnum_chain
from math import isnan

pdbfiles = [ argv[1] ]
prefix = argv[-1]
resnums = []
chains  = []
for arg in argv[2:-1]:
    get_resnum_chain( arg, resnums, chains )

print resnums, chains


def get_pdb_line( lines_out, pdb_lines, resnum_desired, chain_desired ):

    lines = []
    for chain in pdb_lines.keys():

        if ( chain_desired != chain and chain_desired != '' ): continue
        if ( isinstance( resnum_desired, int ) and (resnum_desired not in pdb_lines[ chain ].keys() ) ): continue

        if isinstance( resnum_desired, int ):
            if len( lines ) > 0:
                print 'WARNING! Found residue ', resnum_desired, ' more than once: you may want to specify the chain.'

            for atom_name in pdb_lines[ chain ][ resnum_desired ].keys():
                lines.append( pdb_lines[ chain ][ resnum_desired ][ atom_name ] )

        else: # 'all'
            assert( resnum_desired == 'all' )

            for resnum in pdb_lines[ chain ]:
                for atom_name in pdb_lines[ chain ][ resnum ].keys():
                    lines.append( pdb_lines[ chain ][ resnum ][ atom_name ] )


    if len( lines ) == 0:
        print 'WARNING! Did not find res num ', resnum_desired,
        if chain_desired != '': print ' for chain ', chain_desired,
        print

    for line in lines: lines_out.append( line )


for pdbfile in pdbfiles:
    [ coords, pdb_lines, sequence ] = read_pdb( pdbfile )

    lines_out = []
    for i in range( len( resnums ) ):
        get_pdb_line( lines_out, pdb_lines, resnums[ i ], chains[ i ] )

    pdbfile_out = prefix + pdbfile
    print 'Outputting: ', pdbfile_out
    fid = open( pdbfile_out, 'w' )
    for line in lines_out: fid.write( line+'\n' )
    fid.close()

