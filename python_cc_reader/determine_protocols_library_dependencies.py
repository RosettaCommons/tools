import code_utilities
import library_levels
import pygraph
import inclusion_graph
import sys


if __name__ == "__main__" :

   #inclusion_graph.remove_known_circular_dependencies_from_graph( g )
   #tg = inclusion_graph.transitive_closure( g )

   prot_levels = library_levels.protocols_levels()

   protocols_library_graph = pygraph.digraph()
   projects = code_utilities.rosetta_projects()
   for lib in projects[ "src" ] :
      #print lib
      if lib.find( "protocols" ) != -1 :
         suffix = lib[9:]
         if len( suffix ) == 2 :
            protocols_library_graph.add_node( lib[-1:] ) # e.g. "7" for protocols.7
         elif len( suffix ) == 4 :
            protocols_library_graph.add_node( lib[-1:] + lib[-3:-2] ) # e.g. 4a for protocols_a.4

   #for node in protocols_library_graph.nodes() :
   #    print " protols libray graph node:", node

   #lev, column = library_levels.library_and_column_for_file( prot_levels, "protocols/abinitio/ClassicAbinitio.hh" )
   #print "protocols/abinitio/ClassicAbinitio.hh", str(lev) + str(column)
   #lev, column = library_levels.library_and_column_for_file( prot_levels, "protocols/init/init.hh" )
   #print "protocols/init/init.hh", str(lev) + str(column)

   #sys.exit(0)

   compilable_files, all_includes, file_contents = code_utilities.load_source_tree()
   print "...source tree loaded"
   g = inclusion_graph.create_graph_from_includes( all_includes )
   print "...inclusion graph created"

   for edge in g.edges() :
      file1,file2 = edge
      lib1 = library_levels.libname_for_file( file1 )
      lib2 = library_levels.libname_for_file( file2 )
      if lib1 != "protocols" or lib2 != "protocols" : continue
      lev1,col1 = library_levels.library_and_column_for_file( prot_levels, file1 )
      lev2,col2 = library_levels.library_and_column_for_file( prot_levels, file2 )
      node1 = str(lev1) + col1
      node2 = str(lev2) + col2
      if node2 not in protocols_library_graph.neighbors( node1 ) :
          print "adding edge", node1, node2
          protocols_library_graph.add_edge( node1, node2 )

   protocols_tg = inclusion_graph.transitive_closure( protocols_library_graph )
   for node in protocols_tg.nodes() :
      print "Dependencies for node", node
      for node_neighbor in protocols_tg.neighbors( node ) :
          print "  ", node_neighbor


