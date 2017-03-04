#!/usr/bin/env python

from sys import argv
from os.path import abspath

for file in argv[1:]:
    if ( file.find( '.cc' ) == -1 and file.find( '.hh' ) == -1 ): print 'Skipping because not .cc file: ', file
    lines = open( file ).readlines()
    changed_a_line = False

    new_lines = []
    for line in lines:
        if line.find( 'AUTO-REMOVE' ) > -1:
            changed_a_line = True
        else:
            new_lines.append( line )

    if changed_a_line:
        print 'Changed file: ' , file
        fid = open( file, 'w' )
        for line in new_lines: fid.write( line )
        fid.close()
        continue

