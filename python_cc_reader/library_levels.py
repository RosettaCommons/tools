from python_cc_reader.cpp_parser.code_utilities import *
from python_cc_reader.inclusion_removal.inclusion_graph import *
from optparse import OptionParser
import sys


#Info: Run this script from the source/src directory.  It will check all library levels to make sure there are no illegal dependencies.
   
def initialize_options_parser() :
    parser = OptionParser()
    parser.add_option( "-v", "--verbose", dest="verbose", default=False, action="store_true",
                       help="Include in the list of illegal dependencies those which impact libraries which have not yet been split (i.e. protocols)", )
    return parser

def all_protocols_dirs():
    protdirs = []
    dir = "protocols"
    for item in os.listdir(dir):
        name = os.path.join(dir,item)
        if not os.path.isfile( name ) :
            protdirs.append( item )
            #print "protocols directory:",name
    #print protdirs
    return protdirs


def basic_levels() :
   levels = []
   levels.append( (toplevel_subdirs_of_library( "basic" ),) ) # all the subdirectories of basic are in the same library.
   return levels

def levels_for_library( library_prefix ) :
   projects = rosetta_projects()
   library_sublibs = {}
   sorted_src_projects = sorted( projects[ "src" ] )
   for lib in sorted_src_projects :
       print(lib)
       if lib.find( library_prefix ) != -1 :
           level = int( lib.rpartition(".")[2] )
           print(lib, level)
           if level not in library_sublibs :
               library_sublibs[level] = []
           library_sublibs[level].append( lib )
   levels = []
   nlevels = max( library_sublibs.keys() )
   for level in range(1,nlevels+1) :
      columns = []
      for libname in library_sublibs[ level ] :
         subdirs = toplevel_subdirs_of_library( libname )
         #print libname, subdirs
         #for subdir in subdirs :
         #    print libname, subdir
         columns.append( subdirs )
      levels.append( columns )
   return levels


# take the library levels directly out of the .src.settings files, in case
# someone adds a new top-level directory in core/
def core_levels() :
   return levels_for_library( "core" )

def protocols_levels() :
   return levels_for_library( "protocols" )

def library_and_column_for_file( levels_for_lib, fname ) :
   dirs = fname.split("/")
   if len(dirs) < 3 :
      print("Error in library_and_column_for_file", fname)
      print("Could not find subdirectory of the library directory")
   topleveldir = dirs[1]
   for i in range(len(levels_for_lib)) :
      ilib = levels_for_lib[i]
      for column in range(len(ilib)) :
         if topleveldir in ilib[ column ] :
            if len(ilib) == 1 :
               return i+1, "" # no column name for libraries with only a single entry
            else :
               return i+1, chr( ord('a') + column )
   print("ERROR top level directory",topleveldir,"not found in levels for this library")
   return None, None

class DesiredDependencies :
   def __init__( self ) :
      self.lib_levels_ = [
         ( "platform", [], True ),
         ( "ObjexxFCL", [], True ),
         ( "utility", [], True ),
         ( "numeric", [], True ),
         ( "basic", basic_levels(), True ),
         ( "core", core_levels(), True ),
         ( "protocols", protocols_levels(), True ), # set this to 'True' to look for illegal dependencies in the protocols library
         ( "devel", [], False ),
         ( "apps", [], False ) ]
      for lib in self.lib_levels_ :
          #print lib[0]
          if len( lib[1] ) > 0 :
              count = 0
              for dirs in lib[1] :
                  count += 1
                  #print " ", lib[0] , count, ":",
                  #for directory in dirs :
                  #    print directory,
                  #print
      #print self.lib_levels_[ 6 ][ 0 ]
      #print sum( [ len(x) for x in self.lib_levels_[ 6 ][ 1 ] ] )
      #sys.exit()
   def lib_id( self, libname ) :
      for i in range( self.lib_levels_ ) :
         if libname == self.lib_levels_[ i ][ 0 ]:
            return i
      return -1;

   def treat_illegal_dependency_seriously_for_library( self, libname ) :
      return self.lib_levels_[ libname ][ 2 ]

   def subdir_level( self, libname, subdir_name ) :
      #print "subdir_level", libname, subdir_name
      for i in range( len( self.lib_levels_ )) :
         if libname == self.lib_levels_[ i ][ 0 ] :
            for j in range( len( self.lib_levels_[ i ][ 1 ] )) :
               for k in range( len( self.lib_levels_[ i ][ 1 ][ j ] ) ) :
                  #print "  subdir_level", libname, subdir_name, i, j, k, self.lib_levels_[ i ][ 1 ][ j ][ k ]
                  if subdir_name in self.lib_levels_[ i ][ 1 ][ j ][ k ]:
                     #print "returning", j, k
                     return ( j, k )
      #print "returning -1, -1"
      return ( -1, -1 )

   def levels_for_lib( self, libname ) :
      for i in range( len( self.lib_levels_)):
         if libname == self.lib_levels_[ i ][ 0 ] :
            return self.lib_levels_[ i ][ 1 ]
      return None

   def legal_dependency( self, parent, dependent ) :
      # the dependent is the file doing the #including; the parent is the file being #included
      dep_dirs = dependent.split("/")
      if len( dep_dirs ) < 3 :
         return True
      dep_lib = dep_dirs[ 0 ]
      par_dirs = parent.split("/")
      par_lib = par_dirs[ 0 ]
      dep_libid = self.lib_id( dep_lib )
      par_libid = self.lib_id( par_lib )
      if par_libid == -1 or dep_libid == -1 :
         return True
      if dep_libid < par_libid :
         print("Grossly illegal!", dependent, parent)
         return False
      if dep_libid > par_libid :
         return True
      dep_subdir_level, dep_subdir_column = self.subdir_level( dep_lib, dep_dirs[ 1 ] )
      par_subdir_level, par_subdir_column = self.subdir_level( par_lib, par_dirs[ 1 ] )

      if dep_subdir_id == -1 or par_subdir_id == -1 :
         return False
      if par_subdir_id > dep_subdir_id :
         return False
      if par_subdir_id == dep_subdir_id :
         return par_subdir_column == dep_subdir_column
      return True

   def level_for_lib( self, libname ):
      for i in range( len( self.lib_levels_ )):
         if self.lib_levels_[ i ][0] == libname :
            return i
      return -1

   #def level_for_ns( self, libname, ns ):
   #   for i in xrange( len( self.lib_levels_ )):
   #      if self.lib_levels_[ i ][ 0 ] != libname : continue
   #      for j in xrange( len(self.lib_levels_[ i ][1])):
   #         if ns in self.lib_levels_[ i ][1][ j ]:
   #            return j
   #   return -1


#def purge_illegal_dependencies( dot_lines, liblevels, libname ):
#   newlines = []
#   first_arrow_found = False
#   for line in dot_lines :
#      print line,
#      cols = line.strip().split(" ")
#      if len(cols) != 3 or cols[1] != "->":
#         newlines.append( line )
#      else :
#        if not first_arrow_found :
#          newlines.append( "   compound=true;\n")
#          newlines.append( '   clusterrank="local";\n')
#          for i in xrange(len(liblevels.levels_for_lib( libname))) :
#            newlines.append( "   subgraph cluster_"+libname+str(i)+"{\n" )
#            newlines.append( "      rank = same;\n" )
#            newlines.append( "      label = "+libname+str(i)+";\n" )
#            for ns in levels[i] :
#               newlines.append( '      "' + libname + '::' + ns + '";\n' )
#            newlines.append( "   }\n")
#          first_arrow_found = True
#        d1 = cols[0].split("::")[1].split("\"")[0]
#        d2 = cols[2].split("::")[1].split("\"")[0]
#        l1,col1 = liblevels.subdir_level( libname, d1 )
#        l2,col2 = liblevels.subdir_level( libname, d2 )
#        print d1, l1, d2, l2, col1, col2
#        if l1 == l2 and l1 != -1 and l2 != -1 :
#           #newlines.append( line )
#           newlines.append( "     " + cols[2] + " -> " + cols[0] + '[ dir = "back" ]\n' )
#        if l1 > l2 and l1 != -1 and l2 != -1 :
#           newlines.append( "     " + cols[2] + " -> " + cols[ 0 ] + " [ ltail=cluster_" + libname+str(l2) + ", lhead=cluster_" + libname + str(l1)+ ', dir="back"]\n' )
#   return newlines

def libname_for_file( fname ):
   cols = fname.split("/")
   if len(cols) == 1 :
      return ""
   else:
      return cols[ 0 ]

def subdir_for_file( fname ):
   cols = fname.split("/")
   if len(cols) == 1 or len(cols) == 2 :
      return ""
   return cols[ 1 ]

def find_and_group_illegal_dependencies( ig, count_noncritical_illegals ) :
   # ig = include_graph
   illegal_includes = ( {}, {} ) # element 0 -- highly illegal, element 1 -- regular illegal
   des_deps = DesiredDependencies()
   for edge in ig.edges() :
      n1,n2 = edge
      #print n1, n2
      l1 = libname_for_file( n1 )
      l2 = libname_for_file( n2 )
      if l1 == "" or l2 == "" :
         continue
      llev1 = des_deps.level_for_lib( l1 )
      llev2 = des_deps.level_for_lib( l2 )
      if llev1 > llev2 :
         #perfectly legal
         continue
      elif llev1 < llev2 :
         #if not count_noncritical_illegals and not des_deps.treat_illegal_dependency_seriously_for_library( llev1 ) :
         #   continue
         libdep = (l1,l2)
         if not libdep in illegal_includes[0] : #highly illegal inter-library dependency
            illegal_includes[0][ libdep ] = []
         illegal_includes[0][ libdep ].append( edge )
      elif llev1 == llev2 :
         if not count_noncritical_illegals and not des_deps.treat_illegal_dependency_seriously_for_library( llev1 ) :
            continue

         s1 = subdir_for_file( n1 )
         s2 = subdir_for_file( n2 )
         if s1 == "" or s2 == "" :
            continue
         slev1, col1 = des_deps.subdir_level( l1, s1 )
         slev2, col2 = des_deps.subdir_level( l2, s2 )
         if slev1 > slev2 :
            #perfectly legal
            continue
         elif slev1 == slev2 and col1 == col2:
            #perfectly legal
            continue
         else:
            l1s1 = l1 + "/" + s1
            l2s2 = l2 + "/" + s2
            subdirdep = (l1s1, l2s2)
            if not l1 in illegal_includes[1] :
               illegal_includes[1][l1] = {}
            if not s1 in illegal_includes[1][l1] :
               illegal_includes[1][l1][s1] = {}
            if not subdirdep in illegal_includes[1][l1][s1] :
               illegal_includes[1][l1][s1][ subdirdep ] = []
            illegal_includes[1][l1][s1][ subdirdep ].append( edge )
   return illegal_includes

if __name__ == "__main__" :
   parser = initialize_options_parser()
   (options, args) = parser.parse_args()

   prot_levels = protocols_levels();
   prots_assigned = set([])
   for protset in prot_levels :
       for column in protset:
           #print column
           prots_assigned = prots_assigned.union( column )


   all_prot_levels = set( all_protocols_dirs())
   print(all_prot_levels)
   if all_prot_levels - prots_assigned :
       print("protocol directories not assigned to the hierarchy")
       print(all_prot_levels - prots_assigned)
       print("----------")
   if prots_assigned - all_prot_levels :
       print("protocol directories in the hierarchy that do not exist")
       print(prots_assigned - all_prot_levels)

   #sys.exit()


   #cl = core_levels()
   #cd2_lines = open( "core_d2.dot" ).readlines()
   #cd2_lines2 = purge_illegal_dependencies( cd2_lines, cl, "core" )
   #open( "core_d2_desired2.dot", "w" ).writelines( cd2_lines2 )
   #pl = protocols_levels()
   #pd2_lines = open( "protocols_d2.dot" ).readlines()
   #pd2_lines2 = purge_illegal_dependencies( pd2_lines, pl, "protocols" )
   #open( "protocols_d2_desired.dot", "w" ).writelines( pd2_lines2 )

   compilable_files, all_includes, file_contents = load_source_tree()
   g = create_graph_from_includes( all_includes )
   illegal_includes = find_and_group_illegal_dependencies( g, options.verbose )
   des_deps = DesiredDependencies()
   bad_deps_exist = False
   for key in list(illegal_includes[0].keys()):
      print("Highly illegal dependency class:", key[ 0 ], "dependent on", key[ 1 ])
      for edge in illegal_includes[0][key]:
         print("   " + edge[ 0 ] + " --> " + edge[1])
         bad_deps_exist = True
   for lib in list(illegal_includes[1].keys()):
      for subdir in list(illegal_includes[1][lib].keys()):
         for key in list(illegal_includes[1][lib][subdir].keys()):
            #print key, key[0].partition("/")[2], key[1].partition("/")[2]
            s1,col1 = des_deps.subdir_level( lib, key[0].partition("/")[2])
            s2,col2 = des_deps.subdir_level( lib, key[1].partition("/")[2])
            #print key, s1, s2, col1, col2
            if s1 == -1 or s2 == -1 :
                continue
            if s1 == s2 :
               print("s1 == s2", lib, key[0], key[1], s1, s2, col1, col2)
               assert( col1 != col2 )
               print("Illegal lateral dependency between columns at the same level: directory ", key[0], "dependent on", key[1])
            else:
               print("Illegal intra-library dependency", key[ 0 ], "dependent on", key[ 1 ])
            for edge in illegal_includes[1][lib][subdir][key]:
               print("   " + edge[0] + " --> " + edge[1])
            print()
            bad_deps_exist = True

   if bad_deps_exist :
      sys.exit( 1 )
   else:
      sys.exit( 0 )

