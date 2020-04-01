from python_cc_reader.cpp_parser.code_utilities import *
from python_cc_reader.cpp_parser.code_reader import *

from python_cc_reader.inclusion_removal.inclusion_graph import *
from python_cc_reader.inclusion_removal.test_compile import *
from python_cc_reader.inclusion_removal.inclusion_equivalence_sets import *
from python_cc_reader.inclusion_removal.add_headers import *
from python_cc_reader.inclusion_removal.add_namespaces import *
from python_cc_reader.inclusion_removal.remove_header import *
from python_cc_reader.inclusion_removal.remove_duplicate_headers import *
from python_cc_reader.inclusion_removal.reinterpret_objdump import *
from python_cc_reader.inclusion_removal.dont_remove_include import *

import re
import sys
from python_cc_reader.external.pygraph import pygraph

#from pygraph.algorithms.searching import depth_first_search
import subprocess
import pp
import math


# $1 = num cpu -- default of 8 if none given
# $2 = secret phrase for the parallel-python server; optional

secret_phrase = ""
ncpu = 8

if len(sys.argv) > 1 :
   try :
      ncpu = int(sys.argv[1])
   except :
      print("Could not convert first parameter,", sys.argv[1],"to an integer")
      print("Arguments should be python whole_shebang.py <ncpu> <parallel-python-server-secret>")
      sys.exit(1)
if len(sys.argv) > 2 :
   secret_phrase = sys.argv[2]

# OVERVIEW:
# first create the include graph
# trim this graph to remove known cycles (vector1_bool, vector0_bool)


# second, verify that all files compile on their own

# third, compute the transitive closure graph

# fourth, add all headers in transitive closure

# fifth, compute the equivalence sets, focusing only on files in
# core/ protocols/ devel/ and apps/

# sixth, iterate across all equivalence sets
#   divide up the equivalence sets among the compilable files
#   distribute tasks
#   each node reads the source tree, and computes the transitive closure
#   subgraph, looking only at the files proceeded by the files covered by
#   this task* (*write this subroutine)
#   Then
#   iterating across files in this subset:
#      if the file contains nothing but headers to other headers, continue
#      remove all inessential headers from the file
#      through successive recompiles
#      write out the finished file

# BEGIN CODE

# first create the include graph
# trim this graph to remove known cycles (vector1_bool, vector0_bool)


def try_to_compile_files( list_of_files_to_compile ) :
# second, verify that all files compile on their own
  compilable_files, all_includes, file_contents = load_source_tree()

  any_fail_to_compile = False
  for fname in list_of_files_to_compile :
     if not test_compile_from_lines( expand_includes_for_file( fname, file_contents ) ) :
        print("Error: ", fname, "does not compile on its own")
        any_fail_to_compile = True
  if any_fail_to_compile :
     print("Error: coud not compile all files on their own")
   #sys.exit(0)

if __name__ == "__main__" :

   compilable_files, all_includes, file_contents = load_source_tree()


   ppservers = ()
   job_server = pp.Server(ppservers=ppservers, secret=secret_phrase)


   funcs = ( try_to_compile_files, \
          test_compile_from_lines, central_compile_command, remove_duplicate_headers_from_filelines, no_empty_args, \
          test_compile_extreme, generate_objdump, \
          labeled_instructions, ignore_instructions, regexes_for_instructions, \
          skip_sections, relabel_sections, compare_objdump_lines, \
          add_autoheaders_to_filelines, remove_header_from_filelines, group_and_sort_headers, group_headers, \
          add_using_namespaces_to_lines, remove_using_namespace_from_lines, cleanup_auto_namespace_block_lines, \
          expand_includes_for_file, follow_includes_for_file, \
          auto_ns_comment, transitive_closure, prohibit_removal, trim_unnecessary_headers_from_file, \
          trim_inclusions_from_files_extreme, wrap_compile_extreme, \
          load_source_tree, create_graph_from_includes, remove_known_circular_dependencies_from_graph, \
          known_circular_dependencies, topological_sorting, depth_first_search, total_order_from_graph, \
          transitive_closure, write_file, topologically_sorted, \
          libraries_with_ccfiles_to_examine, scan_compilable_files, libraries_with_hhfiles_to_examine,\
          directories_with_ccfiles_to_examine, directories_with_hhfiles_to_examine, \
          include_for_line, find_all_includes, find_includes_at_global_scope, compiled_cc_files, \
          strip_toendofline_comment )

   modules = ( "re", "subprocess", "code_reader", "pygraph", "subprocess" )

   DRI = dont_remove_include.DontRemoveInclude()

   cf_filtered = list(filter( DRI.attempt_include_removal_for_file, compilable_files ))

   nfiles_to_process = len( cf_filtered )
   nfiles_per_cpu = int( math.ceil( nfiles_to_process / ncpu ) )
   print("Testing compilability with", nfiles_per_cpu, "jobs per cpu")
   sys.stdout.flush()
   cf_subsets = []
   start = 0
   for i in range( ncpu - 1 ) :
      cf_subsets.append( cf_filtered[ start : start + nfiles_per_cpu ] )
      start += nfiles_per_cpu
   cf_subsets.append( cf_filtered[ start : ] )
   jobs = []
   count_jobid = 0
   for cf_subset in cf_subsets :
      count_jobid += 1
      jobs.append( job_server.submit( try_to_compile_files, ( cf_subset, ), funcs, modules ) )
   job_server.wait()
   for job in jobs :
      print(job())


