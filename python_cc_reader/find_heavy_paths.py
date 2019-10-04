from python_cc_reader.cpp_parser.code_utilities import (
    load_source_tree,
    non_fwd_hh_subset,
    cc_subset,
    )
from python_cc_reader.cpp_parser.code_reader import CodeReader
from python_cc_reader.inclusion_removal.inclusion_graph import (
    create_graph_from_includes,
    remove_known_circular_dependencies_from_graph,
    transitive_closure )
from python_cc_reader.inclusion_removal.include_paths import paths_to_destinations
import numpy
import os
import re
import pickle
import blargs


def count_paths_to_ccs_for_headers(source_headers, compilable_files, g):
    sorted_fnames = sorted(compilable_files)
    fname_index = {fname:i for i,fname in enumerate(sorted_fnames)}

    N = len(sorted_fnames)
    edge_matrix = numpy.zeros((N, N), dtype=numpy.int32)
    for i,fname in enumerate(sorted_fnames):
        for neighbor in g.node_neighbors[fname]:
            if neighbor not in fname_index:
                continue
            j = fname_index[neighbor]
            edge_matrix[j,i] = 1

    dest_nodes = numpy.zeros((N,), dtype=numpy.int32)
    for i,fname in enumerate(sorted_fnames):
        if fname[-3:] == ".cc":
            dest_nodes[i] = 1

    n_paths_for_node_total = numpy.zeros((N,), dtype=numpy.int64)
    for source_header in source_headers:
        print("looking at paths for", source_header) 
        print("index for", source_header, fname_index[source_header])
        n_paths_for_node = paths_to_destinations(
            fname_index[source_header], edge_matrix, dest_nodes )
        n_paths_for_node_total += n_paths_for_node.astype(numpy.int64)

    return [(fname, n_paths_for_node_total[i]) for i,fname in enumerate(sorted_fnames)]
        


if __name__ == "__main__":
    with blargs.Parser(locals()) as p:
        p.multiword("focused_headers").cast(lambda x: x.split())
    
    pickle_file = "pickled_heavy_headers.bin"
    if os.path.isfile(pickle_file):
        with open(pickle_file,"rb") as fid:
            graphs = pickle.load(fid)
        compilable_files = graphs["compilable_files"]
        all_includes = graphs["all_includes"]
        file_contents = graphs["file_contents"]
        g = graphs["g"]
    else :
        print("loading source tree")
        compilable_files, all_includes, file_contents = load_source_tree()
        print("creating inclusion graph")
        g = create_graph_from_includes(all_includes)
        remove_known_circular_dependencies_from_graph(g)

        graphs = {}
        graphs["compilable_files"] = compilable_files
        graphs["all_includes"] = all_includes
        graphs["file_contents"] = file_contents
        graphs["g"] = g
        with open(pickle_file,"wb") as fid:
            pickle.dump(graphs, fid)

    for focused_header in focused_headers:
        assert focused_header in file_contents

    header_files = non_fwd_hh_subset(file_contents.keys())
    header_importance = count_paths_to_ccs_for_headers(
        focused_headers, file_contents.keys(), g)
    
    sorted_importance = sorted(header_importance, key=lambda x: x[1])

    for header, importance in sorted_importance:
        print(importance, header)

    
