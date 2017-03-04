import fork_manager
import blargs
import subprocess

from test_compile import *
from code_utilities import *
from inclusion_equivalence_sets import *

import re
import sys

# this class keeps track of which process -- represented by pid -- is responsible for
# which job.  It handles the two callbacks from the ForkManager
class JobManager :
   def __init__( self ) :
      self.jobs = {}
      self.failed_jobs = []
   def success_callback( self, fm, p ) :
      if p not in self.jobs :
         print "Critical error.  Could not find header assigned to process ", p
         for p in self.jobs :
            print "Process ", p, "responsible for", self.jobs[ p ]
         sys.exit(1)
      else :
         del self.jobs[p]
   def error_callback( self, fm, p ) :
      if p not in self.jobs :
         print "Critical error.  Could not find header assigned to process ", p
         for p in self.jobs :
            print "Process ", p, "responsible for", self.jobs[ p ]
         sys.exit(1)
      else :
         self.failed_jobs.append( self.jobs[ p ] )
         del self.jobs[ p ]

def run_job(executable, job):
	 print "Running %s on %s" % (executable, job); sys.stdout.flush()
	 command_list = [ executable, "src/" + job ]
	 return subprocess.call( command_list ) == 0
	
if __name__ == "__main__" :
   with blargs.Parser(locals()) as p :
      p.int( "num_cpu" ).shorthand("n").required()
      p.str( "executable" ).shorthand("e").required()
      
   # includes = { 'test.hh': '' }
   includes = scan_compilable_files()

   re_hh_header  = re.compile("\S*\.(cc|hh)$")
   all_files = includes.keys()
   hh_headers = regex_subset( all_files, re_hh_header )
   jobs = hh_headers

   # print "%d files, running with %d jobs..." % ( len(jobs), num_cpu )

   jm = JobManager()
   fm = fork_manager.ForkManager( num_cpu )
   fm.error_callback = jm.error_callback
   fm.success_callback = jm.success_callback

   for job in jobs:
      pid = fm.mfork()
      if pid == 0 :
         if not run_job(executable, job) : sys.exit(1)
         sys.exit(0)
      else :
         jm.jobs[pid] = job
   fm.wait_for_remaining_jobs()

   print 'Failed jobs:', jm.failed_jobs
