#!/usr/bin/env python

from sys import argv
import string

# git show e9911a757c801ee6be503f068cb9ae00c9d4f46d > general_to_align_diffs.txt
lines = open( 'general_to_align_diffs.txt' ).readlines()
line_num = -1
filename = ''
new_file = False
for line in lines:
    # find all diff lines
    if len( line ) > 4 and line[:4] == 'diff':
        if len( filename ) > 0 and len( patches ) > 0:
            filelines = open( filename ).readlines()
            for patch in patches:
                fileline = filelines[ patch[0]-1 ]
                if ( not fileline == patch[1][1:] ) and ( not fileline == patch[2][1:] ):
                    print
                    print 'PROBLEM'
                    print filename
                    print 'LOOK FOR:'
                    print patch[0],patch[1]
                    print 'ACTUAL:'
                    print fileline

                filelines[patch[0] - 1] = patch[2][1:]
            fid = open( filename, 'w' )
            for fileline in filelines:
                fid.write( fileline )
            fid.close()

        cols = string.split( line )
        filename = cols[3][2:]
        patches = []
        replace_lines = []
        new_file = False

    if len( line ) > 2 and line[:2] == '@@':
        cols = string.split( line )
        line_num = int( string.split( cols[2],',' )[0] ) - 1
        continue

    if len( line ) > 0 and line[0] == '-':
        replace_lines.append( line )
        continue

    if len( line ) > 12 and line[:13] == '--- /dev/null':
        new_file = True
    if new_file: continue

    line_num += 1

    if len( line ) > 0 and line[0] == '+':
        if len( line ) > 2 and line[:3] == '+++': continue
        if ( line.find( 'align' ) > -1 ):
            found_it = False
            for replace_line in replace_lines:
                line_original = replace_line
                line_edit = replace_line.replace( 'general', 'align' )
                line_edit = '+' + line_edit[1:]
                if ( line_edit == line ):
                    found_it = True
                    break
            if not found_it:
                continue
                print '!!!!!!PROBLEM!!!!!!!'
                print filename
                print line
                print replace_lines
            patches.append( [line_num,line,line_original] )


