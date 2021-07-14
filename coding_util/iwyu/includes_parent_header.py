#!/usr/bin/env python

"""Usage: includes_parent_header.py

This script is intended to find instances where a given *.cc or *.hh
doesn't include it's "parent" file as the first include.
That is, if a *.cc doesn't include the identically named *.hh file,
or a *.hh doesn't include the corresponding *.fwd.hh

It will only print results for instance where the parent file actually
exists, so if a corresponding *.hh or *.fwd.hh is omitted, it will be
ignored.

Additionally, the script will complain even if the parent include is present,
but is not the very first header included in the file.

This script is intended to be run from the Rosetta/main/source/src directory,
and will automatically scan all the subdirectories.
"""

from __future__ import print_function

import sys, os, os.path
import glob

skipped_headers = set()

source_directories = ["utility", "numeric", "basic", "core", "protocols", "devel"]
#Skipping ObjexxFCL, because that's really an external library.

def check(fn):
    if fn.endswith('.fwd.hh'):
        return

    if fn.endswith('.cc'):
        parent = fn[:-3] + '.hh'
    elif fn.endswith('.hh'):
        parent = fn[:-3] + '.fwd.hh'
    else:
        print( "UNKNOWN FILETYPE CAPTURED:", fn )
        return

    if not os.path.exists( parent ):
        return

    with open( fn ) as f:
        for line in f:
            if not line.startswith('#include'):
                continue
            if len(line.split()) < 2 or line.split()[1] == '<' + parent + '>':
                return # We found it

    # We went through the entire file and didn't find it.
    print( "Issue with missing include: ", fn )
    print( "#include", '<' + parent + '>' )
    print()

def main():
    # Collect all the *.cc and *.hh files
    for source_dir in source_directories:
        if not os.path.exists( source_dir ):
            print( "ERROR directory", source_dir, "not found - run from the main/source/src/ directory" )
            continue
        for dirpath, dirnames, filenames in os.walk(source_dir):
            for filename in [f for f in filenames if f.endswith(".hh") or f.endswith(".cc")]:
                check( os.path.join(dirpath, filename) )


if __name__ == "__main__":
    if len( sys.argv ) != 1 :
        print( __doc__ )
    else:
        main()
