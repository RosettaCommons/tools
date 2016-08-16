#!/usr/bin/env python

from read_pdb import read_pdb
from sys import argv
from make_tag import make_tag_with_dashes


pdbs = argv[1:]

for main_pdb in pdbs:

    assert( main_pdb[-4:] == '.pdb' )
    outfile = main_pdb.replace( '.pdb', '.REORDER.pdb' )
    fid = open( outfile, 'w' )
    [ coords_main, lines_main, sequence_main ] = read_pdb( main_pdb )

    chains = lines_main.keys()
    chains.sort()
    for chain in chains:
        residues = lines_main[ chain ].keys()
        for residue in residues:
            for atom in lines_main[ chain ][ residue ].keys():
                fid.write(  lines_main[ chain ][ residue ][ atom ]+'\n' )
    fid.close()
    print 'Created: ', outfile
