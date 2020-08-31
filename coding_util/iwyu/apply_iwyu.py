#!/usr/bin/env python

'''apply_iwyu.py - a script to read the output of run_iwyu and apply the listed changed.

Note - this assumes that no changes to the relevant file have happened between the run_iwyu.py
and the current entry. If that's not the case, you could get some odd results.

Run from the Rosetta/main/source directory and provide the script with files or directories,
and the .riwyuf therein will be applied.
'''

from __future__ import print_function

import sys, os
import codecs
import json

from optparse import OptionParser

DEBUG = False

def find_insertion_position(contents):
    # Now we need to find the insertion line.
    # For now, we heuristically put it right after the last preprocessor directive,
    # but before any non-comment/non-preprocessor line.
    insertion_pos = 0
    if_loop_level = 0
    for ii, line in enumerate(contents):
        line = line.strip()
        if len(line) == 0:
            continue
        if line.startswith("// END AUTO HEADERS IWYU"):
            return ii+1 # Use this location
        if line.startswith("//"):
            continue
        if line.startswith("#"):
            if line.startswith('#if') and not line.startswith('#ifndef INCLUDED'):
                if_loop_level += 1
                continue
            if line.startswith('#endif'):
                if_loop_level -= 1
            if if_loop_level == 0:
                insertion_pos = ii+1
        elif if_loop_level == 0:
            break

    return insertion_pos

LIB_LEVELS = {
    "apps":0,
    "devel":1,
    "protocols":2,
    "core":3,
    "basic":4,
    "numeric":5,
    "utility":6,
    "ObjexxFCL":7,
    "cifparse":8,
    "rdkit":9,
    "libxml2":10,
    "libzmq":11,
    "cppdb":12,
    "sqlite3":13,
}
def include_sort_key(entry):
    parts = entry.split('/')
    return ( LIB_LEVELS.get(entry[0], 99), entry )


def apply_changes(filename, instructions, options ):
    with open(filename) as f:
        contents = f.readlines()

    for fn, linenums in instructions['deletions'].items():
        for lineno in linenums:
            lineno = int(lineno) -1  # Compensate for zero-based
            line = contents[lineno]
            splitline = line.split()
            if len(splitline) < 2:
                print("ERROR: Line to remove does not match expected content", filename)
                continue
            if splitline[0] != '#include':
                print("ERROR - Line to remove does not match expected content", filename)
                continue
            if splitline[1][1:-1] != fn:
                print("ERROR:: Line to remove does not match expected content", filename)
                continue
            if "DO NOT AUTO-REMOVE" in line:
                print("Skipping line removal for ", fn, " - comment tells it to stay.", filename)
            contents[lineno] = '// AUTOREMOVED IWYU: ' + line.lstrip()

    insertion_pos = find_insertion_position(contents)

    if insertion_pos == 0 and len(instructions['additions'] != 0 ):
        print("ERROR: Can't find insertion position", filename)
        return # The additions are critical to keep compilability. Don't do modification if we can't add.

    added_lines = [ "\n", "\n", "// AUTO HEADERS IWYU\n" ]
    for fn in sorted( instructions['additions'].keys(), key=include_sort_key ):
        whys = instructions['additions'][fn]
        added_lines.append( '#include <'+fn+"> // AUTO IWYU For " + ' '.join(whys) + '\n' )
    added_lines.append( "// END AUTO HEADERS IWYU\n" )


    contents = contents[:insertion_pos] + added_lines + contents[insertion_pos:]

    with open(filename, 'w') as f:
        f.writelines(contents)

def process_file(filename, options):

    with open(filename) as f:
        instructions = json.load(f)

    apply_changes( filename[:-len('.riwyuf')], instructions, options )

    if not options.nodelete:
        os.remove( filename )

def process_dir(dirname, options):
    for item in os.listdir(dirname):
        name = os.path.join(dirname,item)
        if os.path.isdir(name):
            process_dir(name, options)
        elif os.path.isfile(name):
            if name.endswith(".riwyuf"):
                process_file(name, options)

if __name__ == "__main__":
    if not os.path.exists( "./src" ):
        print( "Script must be run from within Rosetta/main/source/" )
        exit()

    parser = OptionParser(usage="usage: %prog [OPTIONS] [FILES|DIRECTORIES]")
    parser.set_description(__doc__)
    parser.add_option("--nodelete",
      default=False, action="store_true",
      help="Don't delete the .riwyuf file afterwards." )
    parser.add_option("--debug",
      default=False, action="store_true",
      help="Print extra debugging info." )

    (options, args) = parser.parse_args(args=sys.argv[1:])

    DEBUG = options.debug

    for name in args:
        if os.path.isdir(name):
            process_dir(name, options)
        elif os.path.isfile(name) and name.endswith('.riwyuf'):
            process_file(name, options)
        elif os.path.isfile( name+'.riwyuf' ):
            process_file(name+'.riwyuf', options)
        else:
            print("Cannot find file or directory: "+name)
            if os.path.isfile(name):
                print("Corresponding .riwyuf not found")
            exit()
