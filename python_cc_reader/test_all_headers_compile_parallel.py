from python_cc_reader.cpp_parser.code_utilities import *
from python_cc_reader.cpp_parser.code_reader import *
from python_cc_reader.cpp_parser import code_reader

from python_cc_reader.inclusion_removal.inclusion_graph import *
from python_cc_reader.inclusion_removal.test_compile import *
from python_cc_reader.inclusion_removal.inclusion_equivalence_sets import *
from python_cc_reader.inclusion_removal.add_headers import *
from python_cc_reader.inclusion_removal.add_namespaces import *
from python_cc_reader.inclusion_removal.remove_header import *
from python_cc_reader.inclusion_removal.remove_duplicate_headers import *
import re
import sys
from . import pygraph
#from pygraph.algorithms.searching import depth_first_search
import subprocess
import pp
import math

def test_headers_compile( header_list ) :
   error_msgs = []
   for header in header_list :
      if not test_compile( header, devnull=True ) :
         error_msgs.append(  header + " fails to compile on its own" )
   return error_msgs


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

ppservers = ()
try :
   job_server = pp.Server(ppservers=ppservers, secret=secret_phrase)
except:
   print( "Could not connect to parallel-python server.  Is it running?  Did you provide the right secret_phrase?")
   sys.exit(1)

includes = scan_compilable_files()
re_hh_header  = re.compile("\S*\.hh$")
re_hpp_header = re.compile( "\S*\.hpp$")

all_files = list(includes.keys())

hh_headers = regex_subset( all_files, re_hh_header )
hpp_headers = regex_subset( all_files, re_hpp_header )

hh_headers.extend( hpp_headers )
headers = hh_headers


funcs = ( \
                       test_compile, central_compile_command, remove_duplicate_headers_from_filelines, no_empty_args, \
                       test_compile_extreme, generate_objdump, \
                       add_autoheaders_to_filelines, remove_header_from_filelines, group_and_sort_headers, group_headers, \
                       add_using_namespaces_to_lines, remove_using_namespace_from_lines, cleanup_auto_namespace_block_lines, \
                       expand_includes_for_file, follow_includes_for_file, \
                       auto_ns_comment, transitive_closure, prohibit_removal, trim_unnecessary_headers_from_file, \
                       trim_inclusions_from_files_extreme, wrap_compile_extreme, \
                       load_source_tree, create_graph_from_includes, remove_known_circular_dependencies_from_graph, \
                       known_circular_dependencies, topological_sorting, depth_first_search, total_order_from_graph, \
                       transitive_closure, write_file, \
                       libraries_with_ccfiles_to_examine, scan_compilable_files, libraries_with_hhfiles_to_examine,\
                       include_for_line, find_all_includes, find_includes_at_global_scope, compiled_cc_files, \
                       strip_toendofline_comment )

modules = ( "re", "subprocess", "code_reader", "pygraph", "subprocess" )

nfiles_to_process = len( headers )
nfiles_per_cpu = int( math.ceil( nfiles_to_process / ncpu ) )
print("Starting header compilation with", nfiles_per_cpu, "jobs per cpu")
header_subsets = []
start = 0
for i in range( ncpu - 1 ) :
   header_subsets.append( headers[ start : start + nfiles_per_cpu ] )
   start += nfiles_per_cpu
header_subsets.append( headers[ start : ] )
jobs = []
for header_subset in header_subsets :
   jobs.append( job_server.submit( test_headers_compile, ( header_subset, ), funcs, modules ) )

job_server.wait()
output = []
for job in jobs :
   jobout = job()
   if jobout :
      output.extend( jobout )

#print len( output )

if output :
   for outlines in output :
      print(outlines)
   sys.exit(1)
else :
   sys.exit(0)


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

