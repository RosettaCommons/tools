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
import os
import re
import pickle

def count_non_blank_lines(header_name, header_lines):
    re_pound_include = re.compile("#include")

    cr = CodeReader()
    cr.push_new_file(header_name)
    count = 0
    for line in header_lines:
        cr.examine_line(line)
        if re_pound_include.match(cr.commentless_line):
            continue
        if len(cr.commentless_line) == 0:
            continue
        count += 1
    return count

def count_transcluded_nlines(fname, header_nlines, tg):
    return sum(header_nlines[neighb] for neighb in tg.node_neighbors[fname] if neighb in header_nlines) + header_nlines[fname]


if __name__ == "__main__":
    pickle_file = "pickled_heavy_headers.bin"
    if os.path.isfile(pickle_file):
        with open(pickle_file,"rb") as fid:
            graphs = pickle.load(fid)
        compilable_files = graphs["compilable_files"]
        all_includes = graphs["all_includes"]
        file_contents = graphs["file_contents"]
        g = graphs["g"]
        tg = graphs["tg"]
    else :
        print("loading source tree")
        compilable_files, all_includes, file_contents = load_source_tree()
        print("creating inclusion graph")
        g = create_graph_from_includes(all_includes)
        remove_known_circular_dependencies_from_graph(g)
        print("computing transitive closure")
        tg = transitive_closure(g)

        graphs = {}
        graphs["compilable_files"] = compilable_files
        graphs["all_includes"] = all_includes
        graphs["file_contents"] = file_contents
        graphs["g"] = g
        graphs["tg"] = tg
        with open(pickle_file,"wb") as fid:
            pickle.dump(graphs, fid)


    header_files = non_fwd_hh_subset(file_contents.keys())
    print("computing header weights")
    header_nlines = {fname:count_non_blank_lines(fname, file_contents[fname]) for fname in header_files}

    header_expense = {fname:count_transcluded_nlines(fname, header_nlines, tg) for fname in header_files} 

    
    header_importance = [(fname, len(cc_subset(tg.node_incidence[fname])) * header_expense[fname]) for fname in header_files]
    sorted_importance = sorted(header_importance, key=lambda x: x[1])

    for header, importance in sorted_importance:
        print(importance, header)

    
