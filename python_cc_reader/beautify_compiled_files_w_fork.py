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

      includes = scan_compilable_files()
      all_files = includes.keys()

      fm = fork_manager.ForkManager( num_cpu )
      fm.error_callback = hcm.handle_failed_header_compilation
      fm.success_callback = hcm.handle_successful_header_compilation
      for file in all_files :
         pid = fm.mfork()
         if pid == 0 :
            beautifier.beautify_file( file )
            sys.exit(0)

   fm.wait_for_remaining_jobs()
