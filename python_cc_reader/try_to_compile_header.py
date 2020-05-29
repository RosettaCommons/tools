#!/usr/bin/env python

import sys
from .test_compile import test_compile

if len(sys.argv) > 3 or ( not sys.argv[1].endswith('.hh') and not sys.argv[1].endswith('.hpp')) :
   print('Usage: python try_to_compile_header.py <header_name> <optional "skip">')
   print('The header name should either end with .hh or .hpp')
   print('If you are sure that your build settings are correct, an optional third')
   print('argument of "skip" may be given to avoid a quick test compile')
   sys.exit(0)

if len( sys.argv) != 3 or sys.argv[2] != "skip" :
   test_header = "core/types.hh"
   print("Test compiling", test_header, " use 3rd argument 'skip' to avoid this step...")
   sys.stdout.flush()
   if not test_compile( test_header ) :
      test_compile( test_header, True )
      print("Failed to compile test header (", test_header, ").")
      print("Your header may compile, but this script will not report that accurately")
      print("unless you make one or more adjustments.  This script is meant to run within")
      print("the rosetta/rosetta_source/src/ directory.  Are you in that directory currently?")
      print("Have you adjusted the os, nbits, and compiler variables in")
      print("test_compile.py file to reflect your system?")
      print("Check err.log and out.log for messages from your compiler.")
      sys.exit(1)
   else :
      print("...test compilation of", test_header, "passed")

print("Testing compilation of ", sys.argv[1])
sys.stdout.flush()
retc = test_compile( sys.argv[1], True )
if ( retc ) :
   print("Success!")
   sys.exit(0)
else :
   print("Header", sys.argv[1], "failed to compile.  Check err.log and out.log for")
   print("messages from your compiler.")
   sys.exit(1)


