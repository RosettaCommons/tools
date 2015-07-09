#!/usr/bin/python

from read_pdb import read_pdb
from sys import argv
from make_tag import make_tag_with_dashes, make_tag_with_dashes_and_commas


csv = False
if '-csv' in argv:
    argv.pop(argv.index('-csv'))
    csv = True

pdbs = argv[1:]

for main_pdb in pdbs:

    [ coords_main, lines_main, sequence_main, chains, residues ] = read_pdb( main_pdb )

    if csv:
             print make_tag_with_dashes_and_commas( residues )
    else:
             print make_tag_with_dashes( residues )

