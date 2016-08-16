#!/usr/bin/env python

from sys import argv
from os.path import abspath

for file in argv[1:]:
    if ( file.find( '.cc' ) == -1 ): print 'Skipping because not .cc file: ', file
    lines = open( file ).readlines()
    changed_a_line = False

    new_lines = []
    for line in lines:
        if line.find( 'static basic::Tracer' ) > -1 or \
           line.find( 'static THREAD_LOCAL basic::Tracer' ) > -1:

            pos = line.find( 'Tracer' )
            pos_start = pos

            pos = pos+ len( 'Tracer' )
            while line[ pos ] != '"': pos += 1
            pos_open_quote = pos
            pos = pos+1

            num_open_quote = 1
            while pos < len( line ) and line[ pos ] != '"': pos += 1

            pos_close_quote = pos

            tracer_name = line[ pos_open_quote+1:pos_close_quote]
            #print tracer_name

            # compare to what it should be.
            abspathname = abspath( file )
            abspathname = abspathname[ (abspathname.find( 'source/src/' ) + 11 ) : ]
            new_tracer_name = abspathname.replace( '/', '.' ).replace( '.cc','' )
            #print new_tracer_name


            if new_tracer_name != tracer_name:
                changed_a_line = True
                line = line[:pos_open_quote] + '"' + new_tracer_name + '"' + line[pos_close_quote+1:]
                print line

        new_lines.append( line )

    if changed_a_line:
        fid = open( file, 'w' )
        for line in new_lines: fid.write( line )
        fid.close()
        continue

