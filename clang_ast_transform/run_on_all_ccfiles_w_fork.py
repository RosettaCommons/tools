import re
import sys, os
import subprocess

sys.path.insert( 0, os.path.realpath(__file__).rpartition("/")[0]+"/../external" )
sys.path.insert( 0, os.path.realpath(__file__).rpartition("/")[0]+"/../python_cc_reader" )

import blargs

from python_cc_reader.utility import fork_manager
from python_cc_reader.inclusion_removal.test_compile import *
from python_cc_reader.cpp_parser.code_utilities import *
from python_cc_reader.inclusion_removal.inclusion_equivalence_sets import *


# this class keeps track of which process -- represented by pid -- is responsible for
# which job.  It handles the two callbacks from the ForkManager
class JobManager :
   def __init__( self ) :
      self.jobs = {}
      self.failed_jobs = []
   def success_callback( self, fm, p ) :
      if p not in self.jobs :
         print("Critical error.  Could not find header assigned to process ", p)
         for p in self.jobs :
            print("Process ", p, "responsible for", self.jobs[ p ])
         sys.exit(1)
      else :
         del self.jobs[p]
   def error_callback( self, fm, p ) :
      if p not in self.jobs :
         print("Critical error.  Could not find header assigned to process ", p)
         for p in self.jobs :
            print("Process ", p, "responsible for", self.jobs[ p ])
         sys.exit(1)
      else :
         self.failed_jobs.append( self.jobs[ p ] )
         del self.jobs[ p ]

def run_job(executable, job):
	 # print "Running %s on %s" % (executable, job)
	 command_list = executable.split(); command_list.append(  "src/" + job )
	 return subprocess.call( command_list ) == 0
	
if __name__ == "__main__" :
   with blargs.Parser(locals()) as p :
      p.int( "num_cpu" ).shorthand("n").required()
      p.str( "executable" ).shorthand("e").required()
      
   # includes = { 'test.hh': '' }
   os.chdir( "src" )
   jobs = compiled_cc_files()
   os.chdir( ".." )

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

   if jm.jobs:
      print(jm.jobs)
   if jm.failed_jobs:
      print("failed jobs", jm.failed_jobs)

   if jm.jobs or jm.failed_jobs:
      sys.exit(1)
   else:
      sys.exit( 0 )
