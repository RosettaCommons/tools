from python_cc_reader.utility import fork_manager
from python_cc_reader.external.blargs import blargs

from python_cc_reader.cpp_parser.code_utilities import *
from python_cc_reader.inclusion_removal.test_compile import *
from python_cc_reader.inclusion_removal.inclusion_equivalence_sets import *

import re
import sys
#import pygraph
#import subprocess
#import math

# this class keeps track of which process -- represented by pid -- is responsible for
# which header.  It handles the two callbacks from the ForkManager
class HeaderCompilationManager :
   def __init__( self ) :
      self.header_for_job = {}
      self.all_headers_compiled = True
      self.headers_that_failed = []
   def handle_successful_header_compilation( self, fm, p ) :
      if p not in self.header_for_job :
         print("Critical error.  Could not find header assigned to process ", p)
         for p in self.header_for_job :
            print("Process ", p, "responsible for", self.header_for_job[ p ])
         sys.exit(1)
      else :
         del self.header_for_job[p]
   def handle_failed_header_compilation( self, fm, p ) :
      if p not in self.header_for_job :
         print("Critical error.  Could not find header assigned to process ", p)
         for p in self.header_for_job :
            print("Process ", p, "responsible for", self.header_for_job[ p ])
         sys.exit(1)
      else :
         self.headers_that_failed.append( self.header_for_job[ p ] )
         self.all_headers_compiled = False
         del self.header_for_job[ p ]

if __name__ == "__main__" :
   with blargs.Parser(locals()) as p :
      p.int( "num_cpu" ).shorthand("n").default(1)  #required()
      p.multiword('headers').default('').unspecified_default()

   headers = headers.split()

   if not headers:  # headers list was not specified, getting it from source files...
       includes = scan_compilable_files()
       re_hh_header  = re.compile("\S*\.hh$")
       re_hpp_header = re.compile( "\S*\.hpp$")

       all_files = list(includes.keys())

       hh_headers = regex_subset( all_files, re_hh_header )
       hpp_headers = regex_subset( all_files, re_hpp_header )

       hh_headers.extend( hpp_headers )
       headers = hh_headers

   hcm = HeaderCompilationManager()
   fm = fork_manager.ForkManager( num_cpu )
   fm.error_callback = hcm.handle_failed_header_compilation
   fm.success_callback = hcm.handle_successful_header_compilation
   for header in headers :
      pid = fm.mfork()
      if pid == 0 :
         if not test_compile( header, devnull=True ) : sys.exit(1)
         sys.exit(0)
      else :
         #print "assigning header", header, "to process", pid
         hcm.header_for_job[pid] = header
   fm.wait_for_remaining_jobs()

   if hcm.all_headers_compiled :
      sys.exit(0)
   else :
      for h in hcm.headers_that_failed :
         print("Header", h, "fails to compile on its own")
      sys.exit(1)


#all_compile = true
#for header in headers :
#   if not test_compile( header ) :
#      print header, "fails to compile on its own"
#      all_compile = false
#   else :
#      pass
#      #print header, "compiles on its own"

#if all_compile :
#   sys.exit( 0 )
#else :
#   sys.exit( 1 )
