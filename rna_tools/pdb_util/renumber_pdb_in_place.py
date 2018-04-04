#!/usr/bin/env python
from __future__ import print_function
from sys import argv
from os import system
from parse_options import parse_options, get_resnum_chain

def renumber_pdb(pdbnames, new_numbers, chains = [], segids = [], retain_atom_num = 0):

    for pdbname in pdbnames:
        lines = open(pdbname,'r').readlines()

        oldresnum = '   '
        count = 0;
        outid  = open( 'temp.txt','w')
        atomnum  = 0
        newchain = ''
        for line in lines:
            line_edit = line
            if line[0:3] == 'TER':
                continue

            if line_edit[0:4] == 'ATOM' or line_edit[0:6] == 'HETATM':
                #if line[17:20] == 'HOH': continue
                if not (line[16]==' ' or line[16]=='A'): continue
                atomnum += 1
                oldchain = line_edit[21]
                resnum = line_edit[23:26]
                oldsegid = line_edit[72:76]
                if not resnum == oldresnum:
                    count = count + 1
                    oldresnum = resnum

                if ( count <= len( new_numbers ) ):
                    newnum = '%4d' % new_numbers[ count-1 ]
                    if len( chains ) > 0:  newchain = chains[ count - 1]
                    if len( segids ) > 0:  newsegid = segids[ count - 1]
                    if newchain == '': newchain = line_edit[ 22 ]
                    if newsegid == '': newsegid = line_edit[ 72:76 ]
                else:
                    if len( new_numbers) > 0:
                        print('WARNING! residue number %d is greater than length of input numbering %d' % (count, len( new_numbers) ))
                    newnum = '%4d' % count
                    newchain = oldchain
                    newsegid = oldsegid

                if len(newsegid) == 1:
                    newsegid = newsegid+'   '
                elif len(newsegid) == 2:
                    newsegid = newsegid+'  '
                elif len(newsegid) == 3:
                    newsegid = newsegid+' '

                if retain_atom_num:
                    line_edit = '%s%s%s' % (line_edit[0:22], newnum, line_edit[26:] )
                else:
                    #print 'replace: ', line_edit[22:26],'->',newnum
                    line_edit = '%s%5d%s%s%s%s%s%s' % (line_edit[0:6],atomnum,line[11:21],newchain, newnum, line_edit[26:72], newsegid, line_edit[76:] )


                outid.write(line_edit)

        if ( count < len( new_numbers) ): print('WARNING! number of residues %d is less than length of input numbering %d' % (count, len( new_numbers) ))

        outid.close()

        system( 'mv temp.txt '+pdbname )

if __name__ == '__main__':
    retain_atom_num = parse_options( argv, "retain_atom_num", 0 )

    assert( len(argv)>1)

    # Look for integers... perhaps we are specifying particular residue numbers
    pdbnames = []
    new_numbers = []
    chains = []
    segids = []
    for i in range(1, len( argv ) ):
        if argv[i].find( '.pdb' )>0 or not get_resnum_chain( argv[i], new_numbers, chains, segids ):
            pdbnames.append( argv[i] )
    renumber_pdb(pdbnames, new_numbers, chains, segids, retain_atom_num)
