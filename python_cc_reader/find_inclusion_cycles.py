#!/usr/bin/env python

# Copyright (c) 2007-2008 Pedro Matiello <pmatiello@gmail.com>
# License: MIT (see COPYING file)

import sys, re

import pygraph
from pygraph.algorithms.critical import transitive_edges
from pygraph.algorithms.sorting import topological_sorting

splitter = re.compile( "\S+" )

def transitive_closure( g ) :
   rtps = topological_sorting( gr )
   rtps.reverse()
   #print rtps

   reachable = {}
   for n in gr:
      reachable[ n ] = set( gr.node_neighbors[ n ] )

   for n in rtps:
      for neighb in gr.node_neighbors[ n ] :
         reachable[ n ] |= reachable[ neighb ]
         #for r in reachable[ neighb ] :
         #   if not r in reachable[ n ] :
         #      reachable[ n ].add( r )
         for r in reachable[ n ] :
            if not r in gr.node_neighbors[ n ] :
               gr.add_edge( n, r )



def find_cycles(graph):
    """
    Find a cycle in the given graph.

    This function will return a list of nodes which form a cycle in the graph or an empty list if
    no cycle exists.

    @type graph: graph, digraph
    @param graph: Graph.

    @rtype: list
    @return: List of nodes.
    """

    if (type(graph) == pygraph.graph):
        directed = False
    elif (type(graph) == pygraph.digraph):
        directed = True
    else:
        raise InvalidGraphType

    def find_cycle_to_ancestor(node, ancestor):
        """
        Find a cycle containing both node and ancestor.
        """
        path = []
        orignode = node
        while (node != ancestor):
            if (node is None):
                return []
            path.append(node)
            node = spanning_tree[node]
        path.append(node)
        path.reverse()
        return path

    def dfs(node):
        """
        Depht-first search subfunction.
        """
        visited[node] = 1
        # Explore recursively the connected component
        for each in graph[node]:
            if (cycle):
                return
            if (each not in visited):
                spanning_tree[each] = node
                dfs(each)
            else:
                if (directed or spanning_tree[node] is not each):
                    cycle.extend(find_cycle_to_ancestor(node, each))

    visited = {}              # List for marking visited and non-visited nodes
    spanning_tree = {}        # Spanning tree
    cycle = []
    cycles = []

    # Algorithm outer-loop
    for each in graph:
        # Select a non-visited node
        if (each not in visited):
            spanning_tree[each] = None
            root = each
            # Explore node's connected component
            dfs(each)
            if (cycle):
                cycles.append( cycle )
                cycle = []

    return cycles


# Graph creation
gr = pygraph.digraph()
#orig_gr = pygraph.digraph()

# Add nodes and edges

node_filename = sys.argv[1]
edge_filename = sys.argv[2]

node_file = open( node_filename )
lines = node_file.readlines()

for line in lines :
  tok = splitter.findall( line )
  gr.add_node( tok[ 0 ] )
  #orig_gr.add_node( tok[ 0 ] )

node_file.close()
edge_file = open( edge_filename )
lines = edge_file.readlines()

for line in lines :
   toks = splitter.findall( line )
   if len(toks) < 2 :
      print "Ignoring line: ", line
      continue
   gr.add_edge( toks[ 0 ], toks[ 1 ] )
   #orig_gr.add_edge( toks[ 0 ], toks[ 1 ] )

#transitive_closure( gr )

cycles = find_cycles( gr )

if len( cycles ) == 0 :
  print "Found no cycles in inclusion graph"
else :

  count = 0
  for cycle in cycles:
    count += 1
    print "-------"
    print "Cycle #", count
    for node in cycle :
      print node



