from inclusion_graph import *
from code_utilities import *
from code_reader import *
from test_compile import *
from reinterpret_objdump import *
from add_headers import *
from add_namespaces import *
import sys
import os
import subprocess
import pp
import math

# $1 = num cpu
# $2 = secret phrase for the parallel-python server; optional
# $3 .. the cxxtest files that should have headers removed from them; if none given, then all cxxtest files will be used

secret_phrase = ""
ncpu = 8

if len(sys.argv) > 1 :
   try :
      ncpu = int(sys.argv[1])
   except :
      print "Could not convert first parameter,", sys.argv[1],"to an integer"
      print "Arguments should be python whole_shebang.py <ncpu> <parallel-python-server-secret>"
      sys.exit(1)
if len(sys.argv) > 2 :
   secret_phrase = sys.argv[2]
if len(sys.argv) > 3 :
   fnames = sys.argv[3:]
else :
   fnames = compiled_cxxtest_hh_files()


good = True
for fname in fnames :
    if not os.path.isfile( fname ) :
        print "Error, file", fname, "does not exist"
        good = False
if not good :
    sys.exit(1)

nfiles_to_process = len( fnames )

if nfiles_to_process == 1 or ncpu == 1 :
    trim_inclusions_from_cxxtest_files( fnames, 0 )

else:

    ppservers = ()
    try :
        job_server = pp.Server(ppservers=ppservers, secret=secret_phrase)
    except:
        print( "Could not connect to parallel-python server.  Is it running?  Did you provide the right secret_phrase?")
        sys.exit(1)

    # this has got to be refactored
    funcs = ( \
        cxxtest_test_compile, cxxtest_gcc_compile_command, cxxtest_testgen_command, \
        trim_inclusions_from_cxxtest_files, wrap_cxxtest_compile, \
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
    nfiles_per_cpu = int( math.ceil( nfiles_to_process / ncpu ) )

    fnames_subsets = []
    start = 0
    for i in range( ncpu - 1 ) :
        fnames_subsets.append( fnames[ start : start + nfiles_per_cpu ] )
        start += nfiles_per_cpu
    fnames_subsets.append( fnames[ start : ] )
    jobs = []
    count_jobid = 0
    for fnames_subset in fnames_subsets :
        count_jobid += 1
        jobs.append( job_server.submit( trim_inclusions_from_cxxtest_files, ( fnames_subset, count_jobid ), funcs, modules ) )
    job_server.wait()
    for job in jobs :
        print job()

