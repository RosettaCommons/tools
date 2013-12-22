#!/usr/bin/python

from os import popen,system
from os.path import basename,exists
from sys import argv,stderr
import string
from parse_options import parse_options


use_subset = False
subset_residues  = parse_options( argv, "subset", [-1] )
segment_residues = parse_options( argv, "segment", [-1] )
if len( segment_residues ) > 0:
    assert( len( subset_residues ) == 0 )
    assert( 2 * (len(segment_residues)/2) == len(segment_residues ) ) # check even
    for i in range( len(segment_residues)/2):
        for j in range( segment_residues[2*i], segment_residues[2*i+1]+1 ):
            subset_residues.append( j )

excise_residues = parse_options( argv, "excise", [-1] )
if len( excise_residues ) > 0:
    use_excise = True
    pdbfiles = argv[1:-1]
    prefix = argv[-1]
    startseq = 1
    endseq = 10000000000

use_subset = ( len( subset_residues ) > 0 )
use_excise = ( len( excise_residues ) > 0 )

if ( use_excise or use_subset ):
    pdbfiles = argv[1:-1]
    prefix = argv[-1]
    startseq = 1
    endseq = 10000000000
else:
    try:
        pdbfiles = argv[1:-2]
        startseq = int( argv[-2])
        endseq = int( argv[-1])
        prefix = 'region_%d_%d_' % ( startseq, endseq )
    except:
        pdbfiles = argv[1:-3]
        startseq = int( argv[-3])
        endseq = int( argv[-2])
        prefix = argv[-1]

atomnums = []

for pdbfile in pdbfiles:

    if not( exists( pdbfile ) ):
        print "Problem: ",pdbfile, " does not exist!"

    gzipped = 0
    outid = open(prefix+basename(pdbfile),'w')

    if pdbfile[-2:] == 'gz':
        lines = popen('zcat '+pdbfile).readlines()
    else:
        lines = open(pdbfile).readlines()

    i = 0
    oldresidue = '   '
    for line in lines:
        if len( line ) > 4  and line[:4] == 'ATOM':
            currentresidue = line[22:26]
            atomnum = int( line[4:11] )

            if not currentresidue == oldresidue:
                i += 1
            oldresidue = currentresidue

            #if (use_subset and not ( i in subset_residues) ): continue
            try:
                currentresidue_val = int(currentresidue)
                if (use_subset and not ( currentresidue_val in subset_residues) ): continue
                if (use_excise and ( currentresidue_val in excise_residues) ): continue
                if int(currentresidue) < startseq or int(currentresidue) > endseq: continue
            except:
                continue

            atomnums.append( atomnum )

            outid.write(line)

        elif len( line ) > 6 and line[:6] == 'CONECT':
            atomnum1 = int( line[6:11] )
            atomnum2 = int( line[11:16] )
            if (atomnum1 in atomnums) and (atomnum2 in atomnums): outid.write( line )

    outid.close()

    if gzipped:
        command = 'gzip -f '+outfile
        system(command)
