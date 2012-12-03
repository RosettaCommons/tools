#!/usr/bin/env python

from sys import argv
bbdep_orig  = argv[1]

aa1   = { \
        'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F',\
        'GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L',\
        'MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R',\
        'SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}

curr_aa = ''
for line in open( bbdep_orig ) :
    new_aa = 'B3' + aa1[ line[:3] ]
    if curr_aa != new_aa :
        curr_aa = new_aa
        out = open( '%s.rotlib' % new_aa, 'w' )
    out.write( 'UNK' + line[3:] )
