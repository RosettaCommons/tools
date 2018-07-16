from __future__ import print_function
import fork_manager
import blargs
import beautifier
import re
import sys, os

from code_utilities import *

# This script is meant to be run from either the Rosetta/main/source/src/ or 
# the Rosetta/main/source/test/ directories. It reads the scons .settings
# files to determine what files are compiled, and then runs the beautifier
# on all of them.  By default, the script does NOT modify the input files, but
# instead writes the beautified output for X.cc to X.cc.beaut.  Use the
# "overwrite" flag -- this will modify X.cc in place.
#
# Some notes on the beautifier.
#
# 1) Your code will need to compile for the beautifier to work correctly.  If you
# have an un-terminated scope, for instance, the beautifier won't be able to
# parse your source code, and therefore won't be able to tell how deeply any particular
# line should be indented.  That said: while your code must compile, it doesn't need
# to work.
#
# 2) The script is working correctly if the only output messages you see are
# "Preparing to run beautifier on X files" where X is ~10K in the src/ directory
# and quite a lot smaller in the test/ directory.  If the beautifier fails to parse
# a file, or if after beautifying that file, the beautifier determines that the new
# structure of the beautified file is not logically identical to the original file,
# then it will not overwrite the starting file.  It will instead list that file in
# the output as one of the files that it was unable to beautify.  Write me if you
# have created a file that compiles but that the beautifier fails to beautify:
# aleaverfay@gmail.com

# ==================================================================================

# this class keeps track of which process -- represented by pid -- is responsible for
# which file.  It handles the two callbacks from the ForkManager
class FileBeautifierManager :
    def __init__( self ) :
        self.file_for_job = {}
        self.all_files_beautified = True
        self.files_that_failed = []
    def handle_successful_file_beautification( self, fm, pid ) :
        if pid not in self.file_for_job :
            print("Critical error.  Could not find file assigned to process ", pid)
            for pid in self.file_for_job :
               print("Process ", pid, "responsible for", self.file_for_job[ pid ])
            sys.exit(1)
        else :
           del self.file_for_job[pid]
    def handle_failed_file_beautification( self, fm, pid ) :
       if pid not in self.file_for_job :
          print("Critical error.  Could not find file assigned to process ", pid)
          for pid in self.file_for_job :
             print("Process ", pid, "responsible for", self.file_for_job[ pid ])
          sys.exit(1)
       else :
          self.files_that_failed.append( self.file_for_job[ pid ] )
          self.all_files_beautified = False
          del self.file_for_job[ pid ]

# from inclusion_equivalence_sets import *

class Dummy :
    pass

def files_to_beautify() :
    # figure out which directory we're in, either source/src or source/test, and return the list
    # of files to operate on based on that.
    lastdir = os.getcwd().rpartition("/")[2]
    if lastdir == "src" :
        return files_in_src_to_beautify()
    elif lastdir == "test" :
        return files_in_test_to_beautify()
    

def files_in_src_to_beautify() :
    includes = scan_compilable_files()
    all_files = includes.keys()

    all_files = filter( lambda x : x.partition("/")[0] != "ObjexxFCL", all_files )
    all_files.remove( "protocols/noesy_assign/PeakAssignmentOptionKeys.hh" ) # this one doesn't beautify
    return all_files

def files_in_test_to_beautify() :
    all_files = [ x.replace( '../test/', '' ) for x in compiled_cxxtest_hh_files() ]
    return all_files

def beautify_all_files_in_pwd( overwrite, num_cpu, pound_if_setting = "take_if", quiet=False  ) :

    all_files = files_to_beautify()
    return beautify_files_in_parallel( all_files, overwrite, num_cpu, pound_if_setting, quiet )


def beautify_files_in_parallel( file_list, overwrite, num_cpu, pound_if_setting = "take_if", quiet=False ) :

    if not quiet :
        if overwrite :
            print("Preparing to run beautifier on", len(file_list), "files")
        else :
            print("Preparing a dry run of the beautifier on", len(file_list), "files; look for beautified output of X.cc in X.cc.beaut")

    fbm = FileBeautifierManager()
    fm = fork_manager.ForkManager( num_cpu )
    fm.error_callback = fbm.handle_failed_file_beautification
    fm.success_callback = fbm.handle_successful_file_beautification

    opts = beautifier.BeautifierOpts()
    if pound_if_setting :
        opts.pound_if_setting = pound_if_setting

    for fname in file_list :
        pid = fm.mfork()
        if pid == 0 :
            # print("beautifying", fname)
            beautifier.beautify_file( fname, overwrite, opts )
            sys.exit(0)
        else :
            fbm.file_for_job[pid] = fname

    fm.wait_for_remaining_jobs()
    return fbm

def exit_following_beautification( fbm ) :
    if fbm.all_files_beautified :
        sys.exit(0)
    else :
       for fname in fbm.files_that_failed :
           print("File", fname, "could not be beautified")
       sys.exit(1)

if __name__ == "__main__" :
    with blargs.Parser(locals()) as p :
        p.int( "num_cpu" ).shorthand("j").default(1)  #required()
        p.str( "pound_if_setting" ).default("take_if") # can be take_if or take_else
        p.flag("overwrite")

    if not overwrite :
        print("WARNING: In the absence of the --overwrite flag, this script will not change")
        print("any of your source files. It will instead write the beautified output to new")
        print(".beaut files for you to review. Also note: you should not commit these .beaut")
        print("files to git! To actually beautify your source files, run this script with")
        print("the --overwrite flag")
    fbm = beautify_all_files_in_pwd( overwrite, num_cpu, pound_if_setting )
    exit_following_beautification( fbm )
