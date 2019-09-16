import subprocess
import pp
import math
import random
import re
import sys


from python_cc_reader.code_improvement.inclusion_removal.inclusion_graph import *
from python_cc_reader.code_improvement.inclusion_removal.test_compile import *
from python_cc_reader.cpp_parser.code_utilities import *
from python_cc_reader.code_improvement.inclusion_removal.inclusion_equivalence_sets import *
from add_headers import *
from add_namespaces import *
from remove_header import *
from remove_duplicate_headers import *
from code_reader import *
from reinterpret_objdump import *

from python_cc_reader.code_improvement.inclusion_removal.dont_remove_include import *
import python_cc_reader.cpp_parser.code_reader as code_reader
import external.pygraph.pygraph as pygraph

#from pygraph.algorithms.searching import depth_first_search

# $1 = num cpu -- default of 8 if none given
# $2 = secret phrase for the parallel-python server; optional


secret_phrase = ""
ncpu = 8

if __name__ == "__main__":

    if len(sys.argv) > 1 :
       try :
          ncpu = int(sys.argv[1])
       except :
          print("Could not convert first parameter,", sys.argv[1],"to an integer")
          print("Arguments should be python whole_shebang.py <ncpu> <parallel-python-server-secret>")
          sys.exit(1)
    if len(sys.argv) > 2 :
       secret_phrase = sys.argv[2]
    
    ppservers = ()
    try :
       job_server = pp.Server(ppservers=ppservers, secret=secret_phrase)
    except:
       print( "Could not connect to parallel-python server.  Is it running?  Did you provide the right secret_phrase?")
       sys.exit(1)
    
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
    
    compilable_files, all_includes, file_contents = load_source_tree()
    g = create_graph_from_includes( all_includes )
    remove_known_circular_dependencies_from_graph( g )
    
    # second, verify that all files compile on their own
    #any_fail_to_compile = False
    #for fname in compilable_files :
    #   if not test_compile_from_lines( expand_includes_for_file( fname, file_contents ) ) :
    #      print "Error: ", fname, "does not compile on its own"
    #      any_fail_to_compile = True
    #if any_fail_to_compile :
    #   print "Error: coud not compile all files on their own"
    #   sys.exit(0)
    
    # third, compute the transitive closure graph
    tg = transitive_closure( g )
    
    # fourth, add all headers in transitive closure
    tar_everything( "bu_starting_code" )
    add_indirect_headers( tg, compilable_files )
    tar_everything( "bu_transclose_headers_round_0" )
    
    
    # fifth, compute the equivalence sets, focusing only on files in
    # core/ protocols/ devel/ and apps/
    equiv_sets = inclusion_equivalence_sets_from_desired_subgraph( g )
    write_equiv_sets_file( "equivalence_sets_wholeshebang.txt", equiv_sets )
    
    #equiv_sets = read_equiv_sets_file( "equivalence_sets_wholeshebang.txt" )
    
    
    
    funcs = ( \
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
    
    modules = ( "re", "subprocess", "code_reader", "pygraph", "subprocess", "dont_remove_include" )
    
    DRI = dont_remove_include.DontRemoveInclude()
    # sixth, iterate across all equivalence sets
    count_round = 0
    for es in equiv_sets :
       count_round += 1
    
       #TEMP skip rounds 1 through 10
       #if count_round <= 10:
       #   continue
    
       es_filtered = list(filter( DRI.attempt_include_removal_for_file, es ))
       random.shuffle( es_filtered )
       nfiles_to_process = len( es_filtered )
       nfiles_per_cpu = int( math.ceil( nfiles_to_process / ncpu ) )
       print("Starting round", count_round, "with", nfiles_per_cpu, "jobs per cpu")
       sys.stdout.flush()
       es_subsets = []
       start = 0
       for i in range( ncpu - 1 ) :
          es_subsets.append( es_filtered[ start : start + nfiles_per_cpu ] )
          start += nfiles_per_cpu
       es_subsets.append( es_filtered[ start : ] )
       jobs = []
       count_jobid = 0
       for es_subset in es_subsets :
          count_jobid += 1
          jobs.append( job_server.submit( trim_inclusions_from_files_extreme, ( es_subset, count_jobid ), funcs, modules ) )
       job_server.wait()
       for job in jobs :
          print(job())
       #tar_together_files( "bu_wholeshebang_round_" + str( count_round ), compilable_files )
       tar_everything( "bu_wholeshebang_round_" + str( count_round )  )
    
    
