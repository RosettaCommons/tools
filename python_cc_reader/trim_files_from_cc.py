#!/local/bin/python

# The purpose of this script is to traverse a single C.cc file and
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

from . import pygraph
#from pygraph.algorithms.sorting import topological_sorting
from .inclusion_graph import add_using_namespaces, remove_using_namespace, cleanup_auto_namespace_header
from .inclusion_graph import backup_file, restore_backup, cc_subset
from .test_compile import test_compile
from .remove_duplicate_headers import remove_duplicate_headers_from_file
from .add_headers import add_autoheaders_to_file
from .remove_header import remove_header_from_file

def trim_headers_from_cc( C_cc, tg, total_order, id=1 ) :
   print("Examining includes within", C_cc)
   #DEBUG TEMP

   if not test_compile( C_cc, id ) :
      print("Skipping ", C_cc, " since it does not compile as is")
      return

   backup_file( C_cc, "bak" )

   remove_duplicate_headers_from_file( C_cc )

   indirect_headers = []
   for header in tg.node_neighbors[ C_cc ] :
      print(header, tg.edge_label( C_cc, header ), C_cc)
      if tg.edge_label( C_cc, header ) == "indirect" :
         indirect_headers.append( header )

   add_autoheaders_to_file( C_cc, indirect_headers )
   backup_file( C_cc, "tmp" )
   if test_compile( C_cc, id ) :
      print("autoheader addition succeeded")
   else :
      print("Warning: could not resolve compilation issues with ", C_cc)
      restore_backup( C_cc, "bak" )
      return

   backup_file( C_cc, "tmp" )
   namespaces  = add_using_namespaces( C_cc )

   # not all namespaces are visible in all files.  Only add the namespaces
   # that do not break the build.
   if not test_compile( C_cc, id ) :
      restore_backup( C_cc, "tmp" )
      newnamespaces = []
      for ns in namespaces :
         print("Testing the addition of namespace", ns)
         backup_file( C_cc, "tmp" )
         test_namespaces = list( newnamespaces )
         test_namespaces.append( ns )
         print("test namespaces: ")
         for ns2 in test_namespaces :
            print("   ", ns2)
         add_using_namespaces( C_cc, test_namespaces )
         if not test_compile( C_cc, id ) :
            print("could not add namespace ", ns)
            restore_backup( C_cc, "tmp" )
         else :
            newnamespaces.append( ns )
      namespaces = newnamespaces

   for namespace in namespaces :
      print(namespace, " added in auto-namespace block")

   new_headers = []
   #for node in toposort :
   #   if node in tg.node_incidence[ C_cc ] :
   #      new_headers.append( node )

   tot_order_set = []
   for header in tg.node_neighbors[ C_cc ] :
      tot_order_set.append( ( header, total_order[ header ] ) )
   values = sorted( tot_order_set, lambda x, y : cmp( x[1], y[1] ) )
   print("Dependent headers for ", C_cc)
   for val in values :
      print("   ", val[ 0 ], val[ 1 ])
      new_headers.append( val[ 0 ] )

   necessary = [] # tg.node_incidence[ C_cc ]
   for X_hh in new_headers :
      print("...Testing: is ",X_hh, "necessary?", end=' ')
      X_hh_necessary = False
      for deps in tg.node_incidence[ X_hh ] :
         if deps in necessary :
            X_hh_necessary = True
            print("yes.  Dependent of", X_hh, "is necessary.")
            break
      if X_hh_necessary :
         remove_header_from_file( C_cc, X_hh )
         continue
      backup_file( C_cc, "tmp" )
      remove_header_from_file( C_cc, X_hh )
      if not test_compile( C_cc, id ) :
         print("yes.  Removing", X_hh, "breaks the build.")
         restore_backup( C_cc, "tmp" )
         necessary.append( X_hh )
      else :
         print("no.", X_hh, "may be safely removed.")

   all_added_namespaces_removed = True
   for ns in namespaces :
      backup_file( C_cc, "tmp" )
      remove_using_namespace( C_cc, ns )
      print("...Testing: is namespace", ns, "necessary?", end=' ')
      if not test_compile( C_cc, id ) :
         print("yes. using namespace", ns, "must remain")
         restore_backup( C_cc, "tmp" )
         all_added_namespaces_removed = False
      else :
         print("no.  Namespace", ns, " does not need a using declaration") 

   # remove the auto-namespace block if it's empty
   cleanup_auto_namespace_header( C_cc )





