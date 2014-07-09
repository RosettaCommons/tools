#!/usr/bin/python

from read_pdb import read_pdb
from sys import argv
from make_tag import make_tag_with_dashes


pdbs = argv[1:]

for main_pdb in pdbs:

    [ coords_main, lines_main, sequence_main, chains, residues ] = read_pdb( main_pdb )

    print make_tag_with_dashes( residues )

