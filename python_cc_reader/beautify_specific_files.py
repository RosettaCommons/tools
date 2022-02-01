#!/usr/bin/env python

# Run this script from somewhere beneath the source/ directory in the main repository;
# it will look for "source" in the CWD and then move into the root directory for the
# repository to do its business

from beautify_compiled_files_w_fork import *
import subprocess, sys
import argparse
import os.path

if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument("-j","--num_cpu", type=int, default=1)
    parser.add_argument("--pound_if_setting", default="take_if", help="Can be take_if or take_else.")
    parser.add_argument("--dry_run", action="store_true", help="Only attempt to beautify without changing the files.")
    parser.add_argument("-q","--quiet", action="store_true", help="Silent all output produced by this script")
    parser.add_argument("-v","--debug", action="store_true", help="Add additional debugging output.")
    parser.add_argument("file_list", nargs="+", help="The files to beautify.")


    options = parser.parse_args()
    num_cpu = options.num_cpu
    pound_if_setting = options.pound_if_setting
    dry_run = options.dry_run
    quiet = options.quiet
    if options.debug and not options.quiet:
        import python_cc_reader.beauty.beautifier
        python_cc_reader.beauty.beautifier.debug = True

    # First, let's make sure all the files are absolute, so they still work when we change directory
    file_list = [ (fn, os.path.abspath(fn)) for fn in options.file_list ]
    file_list = [ (x, abs) for (x, abs) in file_list if (x[-3:]==".hh" or x[-3:]==".cc") ]
    file_list = [ (x, abs) for (x, abs) in file_list if abs.find("source/src/ui/") == -1 and abs.find( "source/code_templates/" ) == -1  ]


    # ok, let's CD into the root of the git repository
    # which should house the "source" directory
    pwd = os.getcwd();
    pwd_parts = pwd.rpartition( "source" )
    assert( pwd_parts[1] == "source" ) # you must run this script from the main/source directory or one of its subdirectories
    if not quiet :
        print("cd'ing to " + pwd_parts[0])
    os.chdir( pwd_parts[0] )

    if not quiet :
        print("Preparing to beautify: " + ", ".join( [ fn for (fn, abs) in file_list ] ))
    # sys.exit(0)

    fbm = beautify_files_in_parallel( [ abs for (fn, abs) in file_list ], not dry_run, num_cpu, pound_if_setting, quiet )
    exit_following_beautification( fbm )

