#!/local/bin/python

# Argument 1: The inclusion graph
#
# Output: Many .cc files will be modified.  The original file
# C.cc is backed up in C.cc.bak.
#
# The purpose of this script is to traverse all .cc files and
# while removing unnecessary #includes, also
# add back in the minimal set of necessary to compile each
# file. Suppose B.hh #inlcudes A.hh and C.cc #includes B.hh.
# If C.cc depends on the #inclusion of B.hh, it
# could be because B.hh defines a class C.cc uses, or it
# could be because B.hh includes another header A.hh that
# defines a class that C.cc uses.  If C.cc needs A.hh but
# does not need B.hh, then this script will perform "path
# compression" -- it will avoid the #inclusion of B.hh
# and link C.cc directly to A.hh
#
# This script will invoke g++ repeatedly to determine
# which headers C.cc requires.

import pp
from . import pygraph
import _thread
from .pygraph.algorithms.sorting import topological_sorting
from .inclusion_graph import read_inclusion_graph, write_inclusion_graph, transitive_closure
from .inclusion_graph import add_using_namespaces, remove_using_namespace, cleanup_auto_namespace_header
from .inclusion_graph import backup_file, restore_backup, cc_subset, auto_ns_comment
from .test_compile import test_compile
from .remove_duplicate_headers import remove_duplicate_headers_from_file
from .add_headers import add_autoheaders_to_file, group_and_sort_headers
from .remove_header import remove_header_from_file
from .trim_files_from_cc import trim_headers_from_cc

funcs  = ( add_using_namespaces, remove_using_namespace, cleanup_auto_namespace_header, backup_file, restore_backup, cc_subset, test_compile, remove_duplicate_headers_from_file, add_autoheaders_to_file, remove_header_from_file, trim_headers_from_cc, group_and_sort_headers, auto_ns_comment )

modules = ("pygraph", "subprocess", "re" )

def trim_headers_from_ccs( cc_list, tg, total_order, index ) :
   for C_cc in cc_list :
      trim_headers_from_cc( C_cc, tg, total_order, index )

class Submitter :
   def __init__( self, job_server, cc_list ) :
      self.lock = _thread.allocate_lock()
      self.count = 0
      self.job_server = job_server
      self.cc_list = cc_list
      self.jobs = []
      self.batch_size = 50

   def submit_job( self, args = [] ) :
      self.lock.acquire()
      if self.count >= len( self.cc_list ) :
         self.lock.release()
         return
      print("Submitting job ", self.count)
      ccsubset = []
      for i in range( self.batch_size ) :
         if self.count + i == len( self.cc_list ) :
            break
         ccsubset.append( self.cc_list[ self.count + i ] )
      self.jobs.append( self.job_server.submit(trim_headers_from_ccs, ( ccsubset, tg, total_order, self.count), funcs, modules, callback=self.submit_job ) )
      self.count += self.batch_size

      self.lock.release()

def local_pafn( fname ) :
   print(fname)
   return 1123


print("reading inclusion graph...")
g = read_inclusion_graph( "filtered_graph.txt" )
#g = read_inclusion_graph( "medium_graph.txt" )
print("done")
print("computing transitive closure graph...")
tg = transitive_closure( g )
print("done")

for node in tg.nodes() :
   for neighb in tg.node_incidence[ node ] :
      print(neighb, "-->", node)


toposort = topological_sorting( tg )
total_order = {}
count_total_order = 0
for file in toposort :
   count_total_order += 1
   total_order[ file ] = count_total_order
   print(file, count_total_order)

cc_files = cc_subset( tg.nodes() )
# cc_files = group_and_sort_headers( cc_files )

#for C_cc in cc_subset( tg.nodes() ) :
#   trim_headers_from_cc( C_cc, tg, total_order, 1 )


ppservers = ()
job_server = pp.Server(ppservers=ppservers)

submitter = Submitter( job_server, cc_files )

# keep three jobs in the queue...
submitter.submit_job()
submitter.submit_job()
##submitter.submit_job()

#print submitter.jobs[ 1 ].finished

job_server.wait()

for job in submitter.jobs :
  print(job())


