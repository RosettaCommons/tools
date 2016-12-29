#!/usr/bin/env python

from sys import argv
from os.path import abspath

for file in argv[1:]:
    if ( file.find( '.cc' ) == -1 and file.find( '.hh' ) == -1 ): print 'Skipping because not .cc or .hh file: ', file

    abspathname = abspath( file )
    abspathname = abspathname[ (abspathname.find( 'source/src/' ) + 11 ) : ]
    namespace_cols = abspathname.split( '/' )[:-1]

    lines = open( file ).readlines()
    new_lines = []

    in_namespace_block = False
    done_namespace_block = False
    found_and_updated_namespace_block = False

    in_close_namespace_block = False
    done_close_namespace_block = False
    found_and_updated_close_namespace_block = False

    brace_count = 0

    for n in range( len( lines ) ):
        line = lines[ n ]

        if len( line ) > 9 and line[:9] == 'namespace' and line.find( 'AUTO ' ) == -1:
            in_namespace_block = True
        elif in_namespace_block:
            done_namespace_block = True
            in_namespace_block = False
        if done_namespace_block:
            if found_and_updated_namespace_block: # error!
                print 'NAMESPACE BLOCK AGAIN?'
                print file
                print lines[n-1]
                print line
                print new_lines
                exit( 1 )
            for col in namespace_cols:
                new_lines.append( 'namespace %s {\n' % col )
            done_namespace_block = False
            found_and_updated_namespace_block = True

        if brace_count == 0 and len(line)>0 and line[0] == '}':
            in_close_namespace_block = True
        elif in_close_namespace_block:
            done_close_namespace_block = True
            in_close_namespace_block = False

        if done_close_namespace_block and not found_and_updated_close_namespace_block:
            for col in namespace_cols[::-1]:
                new_lines.append( '} //%s\n' % col )
            done_close_namespace_block = False
            found_and_updated_close_namespace_block = True

        #print brace_count,in_close_namespace_block, line[:-1]
        if ( not in_namespace_block and not in_close_namespace_block):
            new_lines.append( line )
            line = line.split( '//' )[0]
            brace_count += line.count( '{' )
            brace_count -= line.count( '}' )
            assert( brace_count >= 0 )

    if not found_and_updated_close_namespace_block:
        for col in namespace_cols[::-1]:
            new_lines.append( '} //%s\n' % col )
        found_and_updated_close_namespace_block = True

    assert( brace_count == 0 )
    if found_and_updated_namespace_block and found_and_updated_close_namespace_block:
        fid = open( file, 'w' )
        for line in new_lines: fid.write( line )
        fid.close()
        continue
    else:
        print 'PROBLEM! ', file, found_and_updated_namespace_block, found_and_updated_close_namespace_block


