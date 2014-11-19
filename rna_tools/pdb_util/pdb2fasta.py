#!/usr/bin/python

from sys import argv,stdout
from os import popen,system
from os.path import exists,basename
from get_sequence import get_sequences
from make_tag import make_tag_with_dashes
removechain = 0
if argv.count('-nochain'):
    removechain = 1

pdbnames = argv[1:]

for pdbname in pdbnames:
    (sequences,chains,resnums) = get_sequences( pdbname, removechain )
    fastaid = stdout
    for i in range( len( sequences ) ):
        fastaid.write( '>%s %s\n' % (basename(pdbname), make_tag_with_dashes(resnums[i], chains[i]) ) )
        fastaid.write( sequences[i] )
        fastaid.write( '\n' )
        fastaid.write( '\n' )
