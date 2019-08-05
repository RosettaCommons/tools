#!/usr/bin/env python

from __future__ import __print_function__
from sys import stdout,argv
from os import system
import string

def extractchain(actualpdbname, out, chains_to_extract):

    if actualpdbname[-3:] =='.gz':
        lines = popen( 'zcat '+actualpdbname).readlines()
    else:
        lines = open(actualpdbname,'r').readlines()

#    out = open(actualpdbname_chain_to_extract,'w')
    for i in range( len(lines)):
        line = lines[i]
        if (line.count('ATOM') or (line.count('HETATM'))  )\
                and (line[21:22] in chains_to_extract ):
            #line = line[0:21]+chain_to_extract+line[22:]
            out.write(line)

actualpdbname = argv[1]
chains_to_extract = argv[2:]

newpdbfile = actualpdbname.replace('.pdb',string.join(chains_to_extract,'')+'.pdb')

with open( newpdbfile, 'w' ) as out:
    print 'Extracting to {}...'.format(newpdbfile)
    extractchain(actualpdbname, out, chains_to_extract)
