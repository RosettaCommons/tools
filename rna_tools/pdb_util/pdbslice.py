#!/usr/bin/env python

from os import popen,system
from os.path import basename,exists
from sys import argv,stderr
import string
from parse_options import parse_options

use_subset = False
[subset_residues,subset_chains, subset_segids] = parse_options( argv, "subset", [[0],['A'],['    ']] )
segment_residues = parse_options( argv, "segments", [-1] )
if len( segment_residues ) > 0:
    assert( len( subset_residues ) == 0 )
    assert( 2 * (len(segment_residues)/2) == len(segment_residues ) ) # check even
    for i in range( int(len(segment_residues)/2) ):
        for j in range( segment_residues[2*i], segment_residues[2*i+1]+1 ):
            subset_residues.append( j )
            subset_chains.append( '' )
            subset_segids.append( '    ' )
[excise_residues, excise_chains, excise_segids] = parse_options( argv, "excise", [[0],['A'],['    ']] )
use_subset = ( len( subset_residues ) > 0 )
use_excise = ( len( excise_residues ) > 0 )
use_range = False
startseq = 0
endseq = 0

if ( use_excise or use_subset ):
    pdbfiles = argv[1:-1]
    prefix = argv[-1]
else:
    use_range = True
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

def matches( res, chain, segid, match_res, match_chain, match_segid ):
    assert( len( match_res ) == len( match_chain ) )
    assert( len( match_res ) == len( match_segid ) )
    for n in range( len( match_res ) ):
        if ( res == match_res[n] and ( match_chain[n] == ''  or match_chain[n] == chain ) and ( match_segid[n] == '    '  or match_segid[n] == segid ) ):
            return True
    return False

for pdbfile in pdbfiles:
    if not( exists( pdbfile ) ):
        print("Problem: ",pdbfile, " does not exist!")

    gzipped = 0
    try:
        outid = open(prefix+basename(pdbfile),'w')
    except:
        print('failed to open outfile')

    if pdbfile[-2:] == 'gz':
        lines = popen('zcat '+pdbfile).readlines()
    else:
        lines = open(pdbfile).readlines()

    for line in lines:
        if len( line ) > 4  and ( line[:6] == 'ATOM  ' or line[:6] == 'HETATM'):
            currentchain = line[21]
            currentresidue = line[22:26]
            currentsegid = line[72:76].strip()
            if len(currentsegid) == 0: currentsegid = '    '
            atomnum = int( line[6:11] )
            try:
                currentresidue_val = int(currentresidue)
                if (use_subset and \
                    not matches( currentresidue_val, currentchain, currentsegid, subset_residues, subset_chains, subset_segids ) ): continue
                if (use_excise and \
                    matches( currentresidue_val, currentchain, currentsegid, excise_residues, excise_chains, excise_segids ) ): continue
                if (use_range and (int(currentresidue) < startseq or int(currentresidue) > endseq) ): continue
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
