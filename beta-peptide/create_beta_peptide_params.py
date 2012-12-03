#!/usr/bin/env python


aa1   = { \
        'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F',\
        'GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L',\
        'MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R',\
        'SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}

aa_list = aa1.keys()

for aa in aa_list:

    new_aa =  'B3' + aa1[aa]

    params_in = aa + '.params'
    params_out = new_aa + '.params'

    lines = open( params_in ).readlines()
    fid = open( params_out, 'w' )

    for line in lines:
        if not (line[:4] == "PROP") : line = line.replace( aa, new_aa )
        if ('AA ' in line): line = line.replace(new_aa, 'UNK')

        if line[:8] == 'ATOM  C ':
            fid.write('ATOM  CM  CAbb CT2  -0.18\n')

        if line[:5]=='LOWER':
            fid.write('ATOM 1HM  Hapo HB   0.09\n')
            fid.write('ATOM 2HM  Hapo HB   0.09\n')

        if line == 'BOND  CA   C  \n':
            fid.write( 'BOND  CA   CM  \n' )
            fid.write( 'BOND  CM   C  \n' )
            fid.write( 'BOND  CM  1HM  \n' )
            fid.write( 'BOND  CM  2HM  \n' )
            continue


        icoor_line = {}
        icoor_line['    N  '] = 'ICOOR_INTERNAL    N      0.000000    0.000000    0.000000   N     CA    CM'
        icoor_line['    CA '] = 'ICOOR_INTERNAL    CA     0.000000  180.000000    1.458001   N     CA    CM\nICOOR_INTERNAL    CM     0.000000   68.800011    1.523258   CA    N     CM'
        icoor_line['    C  '] = 'ICOOR_INTERNAL    C      0.000000   68.800011    1.523258   CM    CA    N'
        icoor_line['  UPPER'] = 'ICOOR_INTERNAL  UPPER  150.000015   63.800018    1.328685   C     CM    CA'
        icoor_line['    O  '] = 'ICOOR_INTERNAL    O   -180.000000   59.199963    1.231016   C     CM  UPPER'
        icoor_line['    CB '] = 'ICOOR_INTERNAL    CB  -122.570366   69.555275    1.529797   CA    N     CM'
        icoor_line['  LOWER'] = 'ICOOR_INTERNAL  LOWER -150.000031   58.300003    1.328685   N     CA    CM\nICOOR_INTERNAL   1HM   119.666321   69.842522    1.089216   CM    CA    C \nICOOR_INTERNAL   2HM   120.751678   69.907356    1.089189   CM    CA   1HM'


        if ( line[:14]=='ICOOR_INTERNAL' ):

            tag = line[14:21]

            if tag in icoor_line.keys():
                fid.write( icoor_line[tag]+'\n' )
                continue


        fid.write( line )

    fid.close()
