#!/usr/bin/env python

from sys import argv


for pdb in argv[1:]:
    lines = open( pdb ).readlines()
    changed_a_line = False
    new_lines = []
    if ( pdb.find( '.py' ) > -1 ): continue
    #print pdb
    for line in lines:
        if line.find( '"VIRT"' ) > -1 and \
           line.find( "atom_type" ) > -1 and \
           line.find( "name" ) > -1 and \
           line.find( "atom_type" ) < line.find( '"VIRT"' ):

            num_changes = 0
            while line.find( '"VIRT"' ) > -1 and\
                  line.find( "atom_type" ) > -1 and \
                  line.find( "name" ) > -1 and \
                  line.find( "atom_type" ) < line.find( '"VIRT"' ) and \
                  num_changes <= 2:
                num_changes += 1
                changed_a_line = True

                pos = line.find( '"VIRT"' )
                pos_virt = pos + len( "VIRT" ) + 1

                # go back two parentheses
                while line[ pos ] != ')': pos -= 1
                pos -= 1
                while line[ pos ] != ')': pos -= 1
                pos_close_paren = pos

                while line[ pos ] != '(': pos -= 1
                pos_open_paren = pos

                atomno = line[ pos_open_paren+1: pos_close_paren ]

                # make our way back to rsd name (behind atom_name)"
                while pos > -1 and line[ pos ] != '.' and line[ pos ] != '>' : pos -= 1
                if pos == -1: print "PROBLEM!",pdb
                pos_period = pos+1

                line = line[:pos_period] + "is_virtual(" + atomno + ") " + line[pos_virt+1:]

            print pdb,':',line[:-1]

        new_lines.append( line )

    if changed_a_line:
        fid = open( pdb, 'w' )
        for line in new_lines: fid.write( line )
        fid.close()
        pass
