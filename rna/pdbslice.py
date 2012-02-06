#!/usr/bin/python

from os import popen,system
from os.path import basename
from sys import argv,stderr
import string



use_subset = 0
subset_residues = []
if argv.count('-subset'):
    use_subset = 1
    pos = argv.index('-subset')
    del argv[pos]

    #stderr.write( 'PDBSLICE using a subset of residues: '  )
    goodint = 1
    while goodint:
        try:
            subset_residue = int(argv[pos])
            subset_residues.append( subset_residue )
            del argv[pos]
            #stderr.write('%d ' % subset_residue )
        except:
            goodint = 0

    #stderr.write( '\n'  )

    pdbfiles = argv[1:-1]

    prefix = argv[-1]
    startseq = 1
    endseq = 10000000000


############################
if argv.count('-segments'):
    use_subset = 1
    pos = argv.index('-segments')
    del argv[pos]

    segment_residues = []
    #stderr.write( 'PDBSLICE using a subset of residues: '  )
    goodint = 1
    while goodint:
        try:
            segment_residue = int(argv[pos])
            segment_residues.append( segment_residue )
            del argv[pos]
            #stderr.write('%d ' % segment_residue )
        except:
            goodint = 0

    #stderr.write( '\n'  )

    pdbfiles = argv[1:-1]

    prefix = argv[-1]
    startseq = 1
    endseq = 10000000000
    for i in range( len(segment_residues)/2):
        for j in range( segment_residues[2*i], segment_residues[2*i+1]+1 ):
            subset_residues.append( j )

    #print subset_residues

use_excise = 0
excise_residues = []
if argv.count('-excise'):
    use_excise = 1
    pos = argv.index('-excise')
    del argv[pos]

    #stderr.write( 'PDBSLICE using a excise of residues: '  )
    goodint = 1
    while goodint:
        try:
            excise_residue = int(argv[pos])
            excise_residues.append( excise_residue )
            del argv[pos]
            #stderr.write('%d ' % excise_residue )
        except:
            goodint = 0

    #stderr.write( '\n'  )

    pdbfiles = argv[1:-1]

    prefix = argv[-1]
    startseq = 1
    endseq = 10000000000


if not use_excise and not use_subset:
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
