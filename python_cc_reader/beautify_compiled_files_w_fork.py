import fork_manager
import blargs
import beautifier
import re
import sys

from code_utilities import *

# from inclusion_equivalence_sets import *

if __name__ == "__main__" :
   with blargs.Parser(locals()) as p :
      p.int( "num_cpu" ).shorthand("n").default(1)  #required()
      p.flag("overwrite")

   includes = scan_compilable_files()
   all_files = includes.keys()
   print "Preparing to run beautifier on", len(all_files), "files"
   
   fm = fork_manager.ForkManager( num_cpu )
   for fname in all_files :
      pid = fm.mfork()
      if pid == 0 :
         print "beautifying", fname
         beautifier.beautify_file( fname, overwrite )
         sys.exit(0)

   fm.wait_for_remaining_jobs()
