# The purpose of this script is to add headers to a set of
# .cc and .hh files that were present in one branch or revision
# that are absent from the current repository.
#
# For example, if ~/rosetta/SVN/mini is the directory where the
# version of mini I'm working with lives, and if
# ~/rosetta/SVN/mini_33042/mini is the directory containing the
# reference copy of mini, and what I want to do is
# inject headers into the working copy that were present
# in the reference copy, then
#
# 1. I would navigate to the ~/rosetta/SVN/mini/ directory
# 2. I would create a list of all the files in src or test
#    that need to be updated with #inclusions that were 
#    present in revision 33042
# 3. I launch this script with two arguments: the name of
#    the list file from step 2, and the directory where
#    revision 33042 lives (~/rosetta/SVN/mini_33042/mini)

import sys, os
from python_cc_reader.cpp_parser.code_utilities import *
from python_cc_reader.inclusion_removal.inclusion_graph import *
from python_cc_reader.inclusion_removal.add_headers import *


class Chdir:         
      def __init__( self, newPath ):  
        self.savedPath = os.getcwd()
        os.chdir(newPath)

      def __del__( self ):
        os.chdir( self.savedPath )

def read_ref_tree( ref_dir ) :
   if ref_dir[ -1 ] == "/" :
      srcdir = ref_dir + "src"
   else :
      srcdir = ref_dir + "/src"
   try :
      newdir = Chdir( srcdir )
   except IOError :
      print("Directory", srcdir, "does not exist!")
      sys.exit(1)
   g = scan_files_to_create_inclusion_graph()
   return g


if len( sys.argv ) != 3 :
   print("Usage: restore_previous_headers.py <files.list> <reference_dir>")
   sys.exit(0)

fnames = open( sys.argv[ 1 ] ).readlines()

ref_g = read_ref_tree( sys.argv[ 2 ] )
ref_tg = transitive_closure( ref_g )

working_g = read_ref_tree( os.getcwd() )
working_tg = transitive_closure( working_g )

print(os.getcwd())

working_total_order = total_order_from_graph( working_g )

for fname in fnames :
   #flines = open( fname.strip() ).readlines()

   try :
      flines = open( fname.strip() ).readlines()
   except IOError :
      print("ERROR: Cannot find file", fname, "Skipping to next file.")
      continue

   inclusions = find_all_includes( fname, flines )
   ref_includes = set( inclusions )
   working_includes = set( inclusions )
   for include in inclusions :
      if include in ref_tg.node_neighbors :
         ref_includes |= set( ref_tg.node_neighbors[ include ] )
      if include in working_tg.node_neighbors :
         working_includes |= set( working_tg.node_neighbors[ include ] )
   to_be_added = []
   missing = []
   working_unreachable = []
   trans_close_included = working_includes
   for include in working_includes :
      if include in working_tg.node_neighbors :
         trans_close_included |= set( working_tg.node_neighbors[ include ] )
   for include in ref_includes :
      if include not in working_includes :
         #print include
         if include in working_g.node_neighbors :
            missing.append( include )
         else :
            working_unreachable.append( include )
   missing = topologically_sorted( working_total_order, missing )
   for miss in missing :
      if miss not in trans_close_included :
         #print miss
         to_be_added.append( miss )
         trans_close_included |= set( working_tg.node_neighbors[ miss ] )
   to_be_added.extend( working_unreachable )
   flines = add_autoheaders_to_filelines( fname, flines, to_be_added )
   #alt_fname = fname.strip() + ".alt"
   #write_file( alt_fname, flines )
   write_file( fname.strip(), flines )

