#!/usr/bin/python

import string
from sys import argv,stderr
from os import popen,system
from os.path import exists,dirname,basename,abspath
from rna_server_conversions import make_rna_rosetta_ready

pdbname = argv[1]

if (pdbname[-4:] != '.pdb' and pdbname[-7:] != '.pdb.gz'):
    pdbname += '.pdb'

removechain = 0
if argv.count('-nochain'):
    pos = argv.index('-nochain')
    del( argv[ pos ] )
    removechain = 1

ignore_chain = 0
if argv.count('-ignorechain'):
    pos = argv.index('-ignorechain')
    del( argv[ pos ] )
    ignore_chain = 1

no_renumber = 0
if argv.count('-no_renumber'):
    pos = argv.index('-no_renumber')
    del( argv[ pos ] )
    no_renumber = 1

chainids = []
if len( argv ) > 2:
    chainids = argv[2:]

if len(chainids) > 0 and len(chainids[0])==1:
    pdbnames = [ pdbname ]
else:
    pdbnames = argv[1:]
    ignore_chain = 1

for pdbname in pdbnames:

    pdblines = string.join( open(pdbname,'r').readlines(), '\n')

    outputstring = make_rna_rosetta_ready( pdblines, removechain, ignore_chain, chainids, no_renumber )

    outfile = string.lower( basename( pdbname ) )
    outfile = outfile.replace( '.pdb', '_RNA.pdb').replace('.gz','');

    outid = open( outfile, 'w')
    outid.write( outputstring )
    print 'Writing ... '+outfile

    outid.close()

