#!/usr/bin/python
import beautify_compiled_files_w_fork
import blargs
import os, sys

# This script is meant to be run from the Rosetta/main/source/ directory.
# It will beautify both the library source files (the src/ directory)
# and the unit test files (the test/ directory).
# This script is working correctly if the only output messages you see say
# "Beautifying src/"
# "preparing to run beautifier on X files",
# "Beautifying test/"
# "preparing to run beautifier on X files".
#
# By default the beautifier does not overwrite the original files, but instead,
# writes out a .beaut file representing the beautified version of your file.
# Use the --overwrite flag to have this script write the beautified output
# on top of the original file.
if __name__ == "__main__" :
    with blargs.Parser(locals()) as p :
        p.int( "num_cpu" ).shorthand("j").default(1)
        p.flag( "overwrite")

    if not overwrite :
        print "WARNING: In the absence of the --overwrite flag, this script will not change"
        print "any of your source files. It will instead write the beautified output to new"
        print ".beaut files for you to review. Also note: you should not commit these .beaut"
        print "files to git! To actually beautify your source files, run this script with"
        print "the --overwrite flag"


    os.chdir( "src" )
    print "Beautifying src/"
    fbm_src = beautify_compiled_files_w_fork.beautify_all_files_in_pwd( overwrite, num_cpu )
    os.chdir( ".." )
    os.chdir( "test" )
    print "Beautifying test/"
    fbm_test = beautify_compiled_files_w_fork.beautify_all_files_in_pwd( overwrite, num_cpu )
    os.chdir( ".." )

    if fbm_src.all_files_beautified and fbm_test.all_files_beautified :
        sys.exit(0)
    else :
        if not fbm_src.all_files_beautified :
            print "Files in src/ that could not be beautified:"
            for fname in fbm_src.files_that_failed :
                print "File", fname, "could not be beautified"
        if not fbm_test.all_files_beautified :
            print "Files in test/ that could not be beautified:"
            for fname in fbm_test.files_that_failed :
                print "File", fname, "could not be beautified"

