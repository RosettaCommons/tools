#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

"""diff5.py Extend a diff3 style git merge conflict, adding a Differ-style delta of the two branches to the output.

The file must have diff3 style conflict marks in it.
Try running `git checkout -m --conflict=diff3 -- <file>` if you have standard 2-way merge conflicts.
"""

from __future__ import print_function

import sys, os
import difflib

class Diffhunk:
    def __init__(self):
        self.leader = []
        self.mine = []
        self.original = []
        self.theirs = []
        self.mymarker = ""
        self.originalmarker = ""
        self.theirmarker = ""
        self.outlines = []

    def parse( self, lines, start ):

        current = self.leader

        pos = start
        while pos < len(lines):
            if lines[pos].startswith("<<<<<<<"):
                self.mymarker = lines[pos]
                current = self.mine
                pos +=1
            elif lines[pos].startswith("|||||||"):
                self.originalmarker = lines[pos]
                current = self.original
                pos += 1
            elif lines[pos].startswith("======="):
                current = self.theirs
                pos += 1
            elif lines[pos].startswith(">>>>>>>"):
                self.theirmarker = lines[pos]
                return pos + 1
            else:
                current.append( lines[pos] )
                pos += 1

        return pos # fallen off the end

    def resolve(self):
        self.outlines.extend( self.leader )
        if len(self.originalmarker) != 0:
            self.outlines.append( self.mymarker )
            self.outlines.extend( self.mine )
            self.outlines.append( "*******\n" )
            self.outlines.extend( difflib.ndiff( self.original, self.mine ) )
            self.outlines.append( self.originalmarker )
            self.outlines.extend( self.original )
            self.outlines.append( "=======\n" )
            self.outlines.extend( difflib.ndiff( self.original, self.theirs ) )
            self.outlines.append( "*******\n" )
            self.outlines.extend( self.theirs )
            self.outlines.append( self.theirmarker )

def main( filename ):
    with open(filename,'r') as f:
        lines = f.readlines()

    pos = 0
    hunks = []
    while pos < len(lines):
        hunk = Diffhunk()
        pos = hunk.parse(lines,pos)
        if len(hunk.mine) != 0 and len(hunk.originalmarker) == 0:
            print( "ERROR: File needs to be in diff3 format" )
            exit()
        hunks.append(hunk)

    with open(filename + '.diff5' ,'w') as f:
        for h in hunks:
            h.resolve()
            f.writelines( h.outlines )

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print( __doc__ )
        exit()
    main( sys.argv[1] )
