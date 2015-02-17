#!/usr/bin/env python

"""Usage: remove_dead_code.py

This script is intended to find and remove dead code in Rosetta.
For this purpose, dead code is considered to be any .cc files which aren't compiled
by the scons system, and any header file which isn't used by a cc file (either
directly or iteratively).

This script is intended to be run from the Rosetta/main/source/ directory.

This script will actually delete the files on disk, so only run from a completely
checked in working directory.

You'll have to check the results manually, as there are some files (PyRosetta, documentation)
which are removed by this script which shouldn't be.
"""

import sys, os
import glob

skipped_headers = set()

source_directories = ["utility", "numeric", "basic", "core", "protocols", "devel"]
#Not apps, because the .src.settings script is funky, and not ObjexxFCL, because that's really an external library.

test_directories = ["utility", "numeric", "basic", "core", "protocols", "devel"]

def collapse_filenames(stem, dictionary):
    names = []
    for key, entry in dictionary.iteritems():
        if type(entry) == list:
            for e in entry:
                names.append( stem + key + "/" + e )
        else:
            print "Cannot interpret: ", entry
    return names

def read_applications():
    app_files = []
    for filename in ("src/apps.src.settings", "src/pilot_apps.src.settings.all"):
        settings = {}
        execfile( filename, settings )
        app_files.extend( collapse_filenames("src/apps/", settings["sources"] ) )
    return [ os.path.normpath(fn + ".cc") for fn in app_files ]

def load_sources(source_dir):
    cc_files = []
    for filename in glob.glob( "src/"+source_dir+"*.src.settings" ):
        settings = {}
        execfile( filename, settings )
        cc_files.extend( collapse_filenames("src/", settings["sources"] ) )
    return [ os.path.normpath(fn + ".cc") for fn in cc_files if not fn.endswith(".cu") ]

def load_tests(test_dir):
    test_files = []
    for filename in glob.glob( "test/"+test_dir+"*.test.settings" ):
        settings = {}
        execfile( filename, settings )
        test_files.extend( collapse_filenames("test/"+test_dir+"/", settings["sources"] ) )
    return [ os.path.normpath(fn + ".cxxtest.hh") for fn in test_files ]

def gather_includes(files):
    includes = set()
    for fn in files:
        if not os.path.exists( fn ):
            print "MISSING: ", fn
            continue
        with open(fn) as f:
            for line in f:
                line = line.strip()
                if not line.startswith("#include"):
                    continue
                if "<" not in line:
                    #print "BAD INCLUDE: ", fn, "->  ", line
                    continue
                name = line[ line.find("<")+1:line.find(">") ]
                if os.path.isfile( "src/" + name ):
                    includes.add( os.path.normpath("src/" + name) )
                elif os.path.isfile( "test/" + name ):
                    includes.add( os.path.normpath("test/" + name) )
                elif os.path.isfile( name ):
                    includes.add( os.path.normpath(name) )
                else:
                    skipped_headers.add( name )
    return includes


def remove_files(directory, keepme):
    #for k in keepme:
    #    if "FuncFactory" in k:
    #        print k
    for dirpath, dirnames, filenames in os.walk( directory ):
        for fn in filenames:
            fullname = dirpath + '/' + fn
            if not ( fn.endswith(".cc") or fn.endswith(".hh") or fn.endswith(".hpp") ):
                #print "IGNORING:", fullname
                continue
            if fullname not in keepme:
                if fullname.endswith(".fwd.hh") and os.path.exists( fullname[:-7] + ".hh" ):
                    #print "FWD MISSING", fullname
                    continue
                if fullname.endswith(".hh") and os.path.exists( fullname[:-3] + ".cc" ):
                    #print "HEADER MISSING", fullname
                    continue
                print "DELETING:", fullname
                os.remove(fullname)

def main():
    cc_files = []
    test_files = []
    header_files = set()
    app_files = read_applications()
    for source_dir in source_directories:
        cc_files.extend( load_sources(source_dir) )
    for test_dir in test_directories:
        test_files.extend( load_tests(test_dir) )

    includes = gather_includes( app_files + cc_files + test_files )
    while includes:
        #print len(includes)
        header_files.update(includes)
        #Find all the includes which we haven't found yet
        includes = gather_includes(includes).difference(header_files)

    for source_dir in source_directories[:-1]: #Don't delete from devel
        remove_files("src/"+source_dir, header_files.union(cc_files + test_files ) )
    #for test_dir in test_directories[:-1]: #Don't delete from devel
    #    remove_files("test/"+test_dir, header_files.union(cc_files + test_files ) )

    #for f in sorted(skipped_headers):
    #    print f

if __name__ == "__main__":
    if len( sys.argv ) != 1 :
        print __doc__
    else:
        main()
