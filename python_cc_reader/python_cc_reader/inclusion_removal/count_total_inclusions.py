from .inclusion_graph import *
from ..cpp_parser.code_utilities import *


def total_inclusion_count():
    g = scan_files_to_create_inclusion_graph()
    remove_known_circular_dependencies_from_graph(g)
    tg = transitive_closure(g)
    count = 0
    for node in tg.nodes():
        count += len(non_fwd_hh_subset(tg.node_neighbors[node]))
        # for non_fwd_hh in non_fwd_hh_subset(tg.node_neighbors[node]):
        #     print(node, "--", tg.edge_label(node, non_fwd_hh), "-->", non_fwd_hh)
    return count

