#!/usr/bin/env python

from os import popen,system
from os.path import dirname,basename,exists
from sys import argv,stderr
import string
from read_pdb import read_pdb
from parse_options import get_resnum_chain
from math import isnan


def Help():
    print argv[0], ' <file.pdb> <desired residues> [prefix] '
    print
    print ' Filenames must end in .pdb.'
    print ' Desired residues can take the form 5-9, C5-9,'
    print '  C:5-9, or C (for all residues in chain c).'
    print " The prefix is prepended on filenames; default is 'subset_'. "
    print
    exit( 0 )

if len( argv ) < 3:
    Help()

pos = 0
resnums = []
chains  = []
segids  = []
pdbfiles = []
prefix = ''
for pos in range( 1, len ( argv) ) :

    # pull out any arguments that look like pdb names
    if argv[ pos ][-4:] == '.pdb':
        pdbfiles.append( argv[pos] )
        continue

    if get_resnum_chain( argv[ pos ], resnums, chains, segids ):
        continue

    if len( prefix ) > 0:
        print 'found two potential prefixes?', prefix, argv[ pos ]
        Help()
    prefix = argv[ pos ]

if len( pdbfiles ) == 0: Help()
if len( resnums ) == 0: Help()
if len( prefix ) == 0: prefix = 'subset_'

print resnums
print chains
print segids


def get_pdb_line( lines_out, pdb_lines, resnum_desired, chain_desired, segid_desired ):

    lines = []
    for chain in pdb_lines.keys():

        if ( chain_desired != chain and chain_desired != '' ): continue
        for segid in pdb_lines[chain].keys():
            if ( segid_desired != segid and segid_desired != '    ' ): continue
            if ( isinstance( resnum_desired, int ) and (resnum_desired not in pdb_lines[ chain ][segid].keys() ) ): continue

            if isinstance( resnum_desired, int ):
                if len( lines ) > 0:
                    print 'WARNING! Found residue ', resnum_desired, ' more than once: you may want to specify the chain.'

                for atom_name in pdb_lines[ chain ][segid][ resnum_desired ].keys():
                    lines.append( pdb_lines[ chain ][segid][ resnum_desired ][ atom_name ] )

            else: # 'all'
                assert( resnum_desired == 'all' )

                for resnum in pdb_lines[ chain ][segid]:
                    for atom_name in pdb_lines[ chain ][segid][ resnum ].keys():
                        lines.append( pdb_lines[ chain ][segid][ resnum ][ atom_name ] )


    if len( lines ) == 0:
        print 'WARNING! Did not find res num ', resnum_desired,
        if chain_desired != '': print ' for chain ', chain_desired,
        print

    for line in lines: lines_out.append( line )


for pdbfile in pdbfiles:
    [ coords, pdb_lines, sequence, _, __, ___] = read_pdb( pdbfile )

    lines_out = []
    for i in range( len( resnums ) ):
        get_pdb_line( lines_out, pdb_lines, resnums[ i ], chains[ i ], segids[ i ] )

    pdbfile_out = '/'.join(filter(None, [
        dirname(pdbfile),
        prefix + basename(pdbfile)
    ]))
    print 'Outputting: ', pdbfile_out
    fid = open( pdbfile_out, 'w' )
    for line in lines_out: fid.write( line+'\n' )
    fid.close()

