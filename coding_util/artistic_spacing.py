#!/usr/bin/env python

from sys import argv
import argparse

#parser = argparse.ArgumentParser(description='Put spaces around equal signs, parentheses, etc. in Rosetta c++ code.')
#parser.add_argument('--force', '-f',help='overwrite files if they previously exist',action="store_true", default=False)
#args = parser.parse_args()
#FORCE = args.force

def Help():
    print argv[ 0 ], ' [ --force ] file1 file2 ...'
    exit( 0 )

if len( argv ) <= 1: Help()

FORCE = False
if '--force' in argv:
    del( argv[ argv.index( '--force' ) ] )
    FORCE = True

operators1_to_space_after   = [ '+', '-', '=', '(', ',', ';','}','{']
operators1_to_space_before  = [ '+', '-', '=', ')' ]

operators1_to_space_after   = [ '+', '-', '=', '(', ',', ';','}','{','<','>']
operators1_to_space_before  = [ '+', '-', '=', ')', '<','>' ]

operators2 = [ '*=', '/=', '==','++','--','+=','-=','!=','<=','>=','()','->','){','};','{}','<<','>>']

operators2_to_space = [ '==','+=','-=','!=','<=','>=']
operators2_do_not_space = [ '--', '++','->', "'-", "-'" ]

ok_spaces = [ ' ',  '\t', '\n' ]
paren_headers = [ 'if','while','for' ]
numchars = [ '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' ]

headers_for_less_than_symbol = [ 'map','vector1','vector','iterator','ptr','Matrix','Vector','cast','pair' ]

def check_headers_for_less_than_symbol( line, pos ):
    for header in headers_for_less_than_symbol:
        if pos > len( header ):
            if line[ pos - len( header ) :  pos  ] == header: return True
    return False

for filename in argv[1:]:

    lines = open( filename ).readlines()
    changed_a_line = False

    new_lines = []
    for line_num in range( len( lines ) ):

        line = lines[ line_num ]

        # fix equal sign
        new_line = ''

        line_length = len( line )

        # don't space in comments or include statements
        max_pos = len( line )
        if line.find( '//' ) > -1: max_pos = line.find( '//' )
        if line.find( '#' ) > -1:  max_pos = min( max_pos, line.find( '#' ) )

        for pos in range( len( line ) ):

            if pos <= max_pos:

                if pos > 0:
                    space_added = False
                    for sign in operators1_to_space_before:
                        if line[ pos ] == sign :
                            if ( line[ pos - 1 ] not in ok_spaces ) and \
                                    line[ pos - 1 : pos + 1 ] not in operators2  and \
                                    line[ pos: pos + 2 ] not in operators2_do_not_space and \
                                    ( sign != '<' or not check_headers_for_less_than_symbol( line, pos ) ) and \
                                     ( sign != '-' or (line[ pos-1 ] != 'e') or (line[ pos-2 ] not in numchars ) ) : # 1e-3

                                new_line += ' '
                                space_added = True
                    if not space_added:
                        for sign in operators2_to_space:
                            if line[ pos : pos + 2] == sign :
                                if ( line[ pos - 1 ] not in ok_spaces ):
                                    new_line += ' '

                if line[ pos ] == '(':
                    for paren_header in paren_headers:
                        if pos >= len( paren_header ) and line[ (pos - len( paren_header )): pos ] == paren_header:
                            new_line += ' '

            new_line += line[ pos ]



            if pos <= max_pos:
                if pos < len( line )-1:
                    space_added = False
                    for sign in operators1_to_space_after:
                        if line[ pos ] == sign:
                            if ( line[ pos + 1 ] not in ok_spaces ) and \
                                    line[ pos : pos + 2 ] not in operators2 and \
                                    line[ pos - 1: pos + 1 ] not in operators2_do_not_space and \
                                    ( sign != '-' or \
                                          ( ( line[ pos-1 ] not in ok_spaces ) and # ( -5 )
                                            ( line[ pos -1 ] not in ['e','E'] or line[ pos - 2] not in numchars ) ) ):
                                    new_line += ' '
                                    space_added = True
                    if not space_added:
                        for sign in operators2_to_space:
                            if line[ pos - 1 : pos + 1 ] == sign :
                                if ( line[ pos + 1 ] not in ok_spaces ):
                                    new_line += ' '

        #new_line = new_line.replace( '}else', '} else' )
        #new_line = new_line.replace( 'else{', 'else {' )
        new_line = new_line.replace( ' ,', ',' )

        # special case:
        new_line = new_line.replace( '" - "', '"-"' )
        new_line = new_line.replace( '( *pose_op )', '(*pose_op)' )

        if new_line != line:
            changed_a_line = True
            print
            print '%s:%d' % ( filename, line_num - 1 )
            print line, new_line,

        new_lines.append( new_line )

    if changed_a_line:
        if FORCE:
            fid = open( filename, 'w' )
            for line in new_lines: fid.write( line )
            fid.close()

if not FORCE:
    print
    print "Changes *not* applied. Rerun with --force flag to apply changes."
else:
    print
    print "Applied changes to files."
