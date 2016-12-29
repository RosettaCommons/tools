#!/usr/bin/env python

from sys import argv


fns = ['Output_boolean','Output_title_text','Output_bool_list','Output_size_list','Output_seq_num_list','output_pair_size','output_pair_size_vector','Output_is_prepend_map']

for pdb in argv[1:]:

    lines = open( pdb ).readlines()
    changed_a_line = False

    new_lines = []
    for line in lines:

        for fn in fns:
            if line.find( fn ) > -1:

                pos = line.find( fn )
                pos_start = pos

                pos = pos+ len( fn )
                while line[ pos ] != '(': pos += 1
                pos_open_paren = pos
                pos = pos+1

                num_open_paren = 1
                while num_open_paren  > 0:
                    if ( line[ pos ] == '(' ): num_open_paren += 1
                    if ( line[ pos ] == ')' ): num_open_paren -= 1
                    #if num_open_paren == 1 and line[ pos ] == ',':
                    #    pos_comma = pos
                    pos += 1

                #assert( pos_comma > 0 )
                pos_close_paren = pos-1

                #print pos_start, pos_open_paren, pos_comma, pos_close_paren

                args = line[ pos_open_paren: pos_close_paren ]
                while args[0] == ' ': args = args[1:]
                while args[-1] == ' ': args = args[:-1]
                if args[-3:] != ' TR' and args[-3:] != 'ing':
                    args += ', TR'
                    changed_a_line = True
                args += ' '

                if args.find( '40, ' ) > -1:
                    args = args.replace(  '40, ', '' )
                    changed_a_line = True
                if args.find( '30, TR' ) > -1:
                    args = args.replace(  '30, TR', 'TR, 30' )
                    changed_a_line = True

                if changed_a_line:
                    line = line[:pos_start] + fn + args + ')' + line[pos_close_paren+1:]
                    print line

        new_lines.append( line )

    if changed_a_line:
        fid = open( pdb, 'w' )
        for line in new_lines: fid.write( line )
        fid.close()
