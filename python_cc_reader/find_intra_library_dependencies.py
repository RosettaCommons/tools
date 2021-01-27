from python_cc_reader.cpp_parser import code_utilities
from python_cc_reader.inclusion_removal import library_levels
from python_cc_reader.inclusion_removal import pygraph
from python_cc_reader.inclusion_removal import inclusion_graph
from python_cc_reader.inclusion_removal import determine_protocols_library_dependencies
import sys


# $1 == the name of the src.settings file to be read.

if __name__ == "__main__" :

   srcsettings_fname = sys.argv[1]
   f = open( srcsettings_fname ).read()
   exec( f )
   topleveldir_graph = pygraph.digraph()
   libname = ""
   for dirname in sources :
      cols = dirname.split("/")
      if libname == "" :
         libname = cols[0]
      else:
         assert( libname == cols[0] )
      if len(cols) < 2 : continue
      if cols[1] not in topleveldir_graph.nodes() :
         topleveldir_graph.add_node( cols[1] )

   print(libname, srcsettings_fname)
   liblevel, libcolumn = determine_protocols_library_dependencies.library_level_and_column_from_src_settingsfile_name( libname, srcsettings_fname )

   prot_levels = library_levels.protocols_levels()

   compilable_files, all_includes, file_contents = code_utilities.load_source_tree()
   print("...source tree loaded")
   g = inclusion_graph.create_graph_from_includes( all_includes )
   print("...inclusion graph created")

   topleveldir_graph_nodes = topleveldir_graph.nodes() # save this for rapid lookup

   for edge in g.edges() :
      file1,file2 = edge
      lib1 = library_levels.libname_for_file( file1 )
      lib2 = library_levels.libname_for_file( file2 )
      if lib1 != libname or lib2 != libname : continue

      lev1,col1 = library_levels.library_and_column_for_file( prot_levels, file1 )
      if lev1 != liblevel or col1 != libcolumn : continue
      lev2,col2 = library_levels.library_and_column_for_file( prot_levels, file2 )
      if lev2 != liblevel or col2 != libcolumn : continue
      node1 = file1.split("/")[1]
      node2 = file2.split("/")[1]
      if node1 not in topleveldir_graph_nodes:
         print("node1:", node1)
         print("all nodes:", topleveldir_graph_nodes)
         assert( node1 in topleveldir_graph_nodes )
      if node2 not in topleveldir_graph_nodes:
         print("node2:", node2)
         print("all nodes:", topleveldir_graph_nodes)
         assert( node2 in topleveldir_graph_nodes )

      if node2 not in topleveldir_graph.neighbors( node1 ) :
          print("adding edge", node1, node2)
          topleveldir_graph.add_edge( node1, node2 )

   #_tg = inclusion_graph.transitive_closure( protocols_library_graph )
   for node in topleveldir_graph.nodes() :
      print("Dependencies for node", node)
      for node_neighbor in topleveldir_graph.neighbors( node ) :
          print("  ", node_neighbor)


