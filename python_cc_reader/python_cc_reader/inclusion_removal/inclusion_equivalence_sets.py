
import sys
from ..external.pygraph import pygraph

def inclusion_equivalence_sets(inclusion_graph):
    dg = pygraph.digraph()
    dg.add_nodes(inclusion_graph.nodes())
    for edge in inclusion_graph.edges():
        dg.add_edge(edge[0], edge[1])

    for node in dg.nodes():
        for neighb in dg.node_neighbors[node]:
            if node == neighb:
                print("Node connects to itself?!:", node)
                sys.exit(1)

    start = "Start"
    color = 0
    red = "red"
    blue = "blue"

    dg.add_node(start)
    remaining = []
    for node in dg.nodes():
        remaining.append(node)
        if node != start and len(dg.node_neighbors[node]) == 0:
            # print( "connecting node ", node, "to the start node")
            dg.add_edge(node, start)
        dg.add_node_attribute(node, red)

    active_nodes = []  # the set of blue nodes with red children
    node_equivalence_sets = []  # the return value, a list of lists.

    dg.node_attr[start][color] = blue
    active_nodes.append(start)
    last_len_remaining = len(remaining)
    remaining.remove(start)

    count = 0
    while len(remaining) != 0:
        count += 1
        len_remaining = len(remaining)
        if last_len_remaining == len_remaining:
            print("Critical error: colored no nodes blue in the previous round!")
            for remnant in remaining:
                print(remnant)
            break
        print(count, len_remaining, last_len_remaining)
        last_len_remaining = len_remaining
        tried_this_round = []
        newly_blued = []
        newly_inactive = []
        for active in active_nodes:
            # print " ", active
            actives_children_all_blue = True
            for node_incident in dg.node_incidence[active]:
                if dg.node_attr[node_incident][color] == blue:
                    continue
                if node_incident in tried_this_round:
                    continue
                # print "  ", node_incident
                tried_this_round.append(node_incident)
                turn_this_node_blue = True
                for parent in dg.node_neighbors[node_incident]:
                    # print "   ", parent, dg.node_attr[ parent ][ color ]
                    if dg.node_attr[parent][color] == red:
                        turn_this_node_blue = False
                        actives_children_all_blue = False
                        break
                if turn_this_node_blue:
                    # print("turning node", node_incident, "blue in round", count)
                    newly_blued.append(node_incident)
            if actives_children_all_blue:
                newly_inactive.append(active)
        for inactive in newly_inactive:
            active_nodes.remove(inactive)
        for new_blue in newly_blued:
            # print count, new_blue
            dg.node_attr[new_blue][color] = blue
            remaining.remove(new_blue)
            active_nodes.append(new_blue)
        if len(newly_blued) == 0:
            print("Error:: no new blue nodes found!")
            for blue_node in active_nodes:
                print("blue node:", blue_node)
                for node_incident in dg.node_incidence[active]:
                    print("  child", node_indicent)
                    for parent in dg.node_neighbors[node_incident]:
                        if dg.node_attr[parent][color] == red:
                            print("   red parent:", parent)
                        
            for remnant in remaining:
                print("remaining node:", remnant)
            break
        node_equivalence_sets.append(newly_blued)

    return node_equivalence_sets
