#!/usr/bin/python

from read_pdb import read_pdb
from sys import argv
from make_tag import make_tag_with_dashes


pdbs = argv[1:]

for main_pdb in pdbs:

    [ coords_main, lines_main, sequence_main ] = read_pdb( main_pdb )

    chains = lines_main.keys()
    chains.sort()
    for chain in chains:
        residues = lines_main[ chain ].keys()
        residues.sort()
        rsd_num = []
        for m in residues: rsd_num.append( m )
        print make_tag_with_dashes( rsd_num )
