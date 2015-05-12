import fork_manager
import blargs
import beautifier
import re
import sys

from code_utilities import *

# this class keeps track of which process -- represented by pid -- is responsible for
# which file.  It handles the two callbacks from the ForkManager
class FileBeautifierManager :
    def __init__( self ) :
        self.file_for_job = {}
        self.all_files_beautified = True
        self.files_that_failed = []
    def handle_successful_file_beautification( self, fm, pid ) :
        if pid not in self.file_for_job :
            print "Critical error.  Could not find file assigned to process ", pid
            for pid in self.file_for_job :
               print "Process ", pid, "responsible for", self.file_for_job[ pid ]
            sys.exit(1)
        else :
           del self.file_for_job[pid]
    def handle_failed_file_beautification( self, fm, pid ) :
       if pid not in self.file_for_job :
          print "Critical error.  Could not find file assigned to process ", pid
          for pid in self.file_for_job :
             print "Process ", pid, "responsible for", self.file_for_job[ pid ]
          sys.exit(1)
       else :
          self.files_that_failed.append( self.file_for_job[ pid ] )
          self.all_files_beautified = False
          del self.file_for_job[ pid ]

# from inclusion_equivalence_sets import *

class Dummy :
    pass

if __name__ == "__main__" :
    with blargs.Parser(locals()) as p :
        p.int( "num_cpu" ).shorthand("n").default(1)  #required()
        p.str( "pound_if_setting" ).default("take_if") # can be take_if or take_else
        p.flag("overwrite")

    includes = scan_compilable_files()
    all_files = includes.keys()
    print "Preparing to run beautifier on", len(all_files), "files"

    fbm = FileBeautifierManager()
    fm = fork_manager.ForkManager( num_cpu )
    fm.error_callback = fbm.handle_failed_file_beautification
    fm.success_callback = fbm.handle_successful_file_beautification

    opts = beautifier.BeautifierOpts()
    if pound_if_setting :
        opts.pound_if_setting = pound_if_setting

    for fname in all_files :
        pid = fm.mfork()
        if pid == 0 :
            print "beautifying", fname
            beautifier.beautify_file( fname, overwrite, opts )
            sys.exit(0)
        else :
            fbm.file_for_job[pid] = fname

    fm.wait_for_remaining_jobs()
    if fbm.all_files_beautified :
        sys.exit(0)
    else :
       for fname in fbm.files_that_failed :
           print "File", fname, "could not be beautified"
       sys.exit(1)
