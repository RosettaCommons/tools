from python_cc_reader.inclusion_removal.inclusion_graph import *
from python_cc_reader.cpp_parser.code_utilities import *

# import os

if __name__ == "__main__":
    cxxtest_files = compiled_cxxtest_hh_files()
    for cxxtest_file in cxxtest_files:
        print("going to add transitive includes to ", cxxtest_file)

    compilable_files, all_includes, file_contents = load_source_tree()
    g = create_graph_from_includes(all_includes)
    remove_known_circular_dependencies_from_graph(g)

    tg = transitive_closure(g)

    # load the set of cxxtest files that need to have headers
    # included
    for cxxtest_file in cxxtest_files:
        add_indirect_headers_to_file_outside_tg(tg, cxxtest_file)
