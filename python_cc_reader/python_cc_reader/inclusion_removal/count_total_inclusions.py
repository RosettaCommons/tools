from .inclusion_graph import *
from ..cpp_parser.code_utilities import *


def total_inclusion_count():
    g = scan_files_to_create_inclusion_graph()
    remove_known_circular_dependencies_from_graph(g)
    tg = transitive_closure(g)
    count_all_nonfwd = 0
    count_internal_nonfwd = 0
    for node in cc_subset(tg.nodes()):
        count_internal_nonfwd += len(non_fwd_hh_subset(tg.node_neighbors[node]))
        count_all_nonfwd += len(regex_exclude(tg.node_neighbors[node], fwdhh_regex()))
        # for non_fwd_hh in regex_exclude(tg.node_neighbors[node], fwdhh_regex()):
        #     print(node, "--", tg.edge_label(node, non_fwd_hh), "-->", non_fwd_hh)
    return count_internal_nonfwd, count_all_nonfwd

