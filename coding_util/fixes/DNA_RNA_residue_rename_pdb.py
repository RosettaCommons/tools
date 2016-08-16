#!/usr/bin/env python

from sys import argv


for pdb in argv[1:]:
    assert( pdb[-4:] == '.pdb' )

    lines = open( pdb ).readlines()

    resmap = { '   A':'  DA', '   C':'  DC','   T':'  DT','   G':'  DG',
               '  rA':'   A', '  rC':'   C','  rU':'   U','  rG':'   G' }
    new_lines = []
    for line in lines:
        if len( line ) > 20 and line[:4] == 'ATOM':
            resname = line[16:20]
            if resname in resmap.keys():
                line = line[:16] + resmap[ resname ] + line[20:]
        new_lines.append( line )

    fid = open( pdb, 'w' )
    for line in new_lines: fid.write( line )
    fid.close()
