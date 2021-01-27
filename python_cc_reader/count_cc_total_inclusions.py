from pythonc_cc_reader.inclusion_removal.inclusion_graph import *
from pythonc_cc_reader.cpp_parser.code_utilities import *


def total_inclusion_count():
    g = scan_files_to_create_inclusion_graph()
    tg = transitive_closure(g)
    count = 0
    for node in cc_subset(tg.nodes()):
        count += len(non_fwd_hh_subset(tg.node_neighbors[node]))
        for non_fwd_hh in non_fwd_hh_subset(tg.node_neighbors[node]):
            print(node, " ", non_fwd_hh)
    return count


if __name__ == "__main__":
    print(total_inclusion_count())
