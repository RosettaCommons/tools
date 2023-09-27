#!/usr/bin/env python

from read_pdb import read_pdb
from sys import argv
from make_tag import make_tag_with_dashes


pdbs = argv[1:]

for main_pdb in pdbs:

    assert( main_pdb[-4:] == '.pdb' )
    outfile = main_pdb.replace( '.pdb', '.REORDER.pdb' )
    fid = open( outfile, 'w' )
    [ coords_main, lines_main, sequence_main, chains_main, residues_main, segids_main ] = read_pdb( main_pdb )

    chains = list(lines_main.keys())
    chains.sort()
    for chain in chains:
        segids = list(lines_main[ chain ].keys())
        segids.sort()
        for segid in segids:
            residues = list(lines_main[ chain ][ segid ].keys())
            residues.sort()
            for residue in residues:
                for atom in list(lines_main[ chain ][ segid ][ residue ].keys()):
                    fid.write(  lines_main[ chain ][ segid ][ residue ][ atom ]+'\n' )
    fid.close()
    print('Created: ', outfile)
