#!/usr/bin/env python

from sys import argv


for pdb in argv[1:]:

    lines = open( pdb ).readlines()
    changed_a_line = False

    new_lines = []
    for line in lines:
        while line.find( 'Contain_seq_num' ) > -1:

            changed_a_line = True

            pos = line.find( 'Contain_seq_num' )
            pos_start = pos

            pos = pos+ len( 'Contain_seq_num' )
            while line[ pos ] != '(': pos += 1
            pos_open_paren = pos
            pos = pos+1

            num_open_paren = 1
            pos_comma = 0
            while num_open_paren  > 0:
                if ( line[ pos ] == '(' ): num_open_paren += 1
                if ( line[ pos ] == ')' ): num_open_paren -= 1
                if num_open_paren == 1 and line[ pos ] == ',':
                    assert( pos_comma == 0 )
                    pos_comma = pos
                pos += 1

            pos_close_paren = pos-1

            #print pos_start, pos_open_paren, pos_comma, pos_close_paren

            list_name = line[ pos_comma+1:pos_close_paren]
            while list_name[0] == ' ': list_name = list_name[1:]
            while list_name[-1] == ' ': list_name = list_name[:-1]
            if list_name[-2:] == '()': list_name = ' ('+list_name+')'

            line = line[:pos_start] + list_name + '.has_value(' + line[ pos_open_paren+1:pos_comma] + ')' + line[pos_close_paren+1:]
            print line

        new_lines.append( line )

    if changed_a_line:
        fid = open( pdb, 'w' )
        for line in new_lines: fid.write( line )
        fid.close()
