import sys

from .remove_headers_from_all_files import (
    create_inclusion_graph_from_source,
    add_transitive_closure_headers,
    remove_includes_from_files_in_parallel
)
from .inclusion_graph import (
    create_graph_from_includes,
    remove_known_circular_dependencies_from_graph,
    transitive_closure,
    add_indirect_headers,
    inclusion_equivalence_sets_from_desired_subgraph,
    write_equiv_sets_file,
    trim_unnecessary_headers_from_file,
    total_order_from_graph,
    wrap_compile_extreme,
    wrap_compile_w_surrogates,
    prohibit_removal_of_non_subgraph_headers
)
from .dont_remove_include import DontRemoveInclude

from .test_compile import test_compile, tar_everything, generate_objdump_for_file
from .add_headers import add_autoheaders_to_file, write_file
from ..external.pygraph import pygraph


def add_inclusions_for_subgraph(tg, starting_files, file_contents):
    for fname in tg.node_neighbors:
        transcludes = set([])
        for starting_file in starting_files:
            if starting_file in tg.node_neighbors[fname]:
                if tg.edge_label(fname, starting_file) == "indirect":
                    transcludes.add(starting_file)
                for transclude in tg.node_neighbors[starting_file]:
                    if tg.edge_label(fname, transclude) == "indirect":
                        transcludes.add(transclude)
        if len(transcludes) > 0:
            # print("adding autoheaders to", fname)
            add_autoheaders_to_file(fname, list(transcludes))
            file_contents[fname] = open(fname).readlines()

def subgraph_from_starting_nodes(g, tg, starting_files):
    subgraph = pygraph.digraph()
    for fname in starting_files:
        subgraph.add_node(fname)
        # print("adding subgraph node", fname)
        for descendent in tg.node_incidence[fname]:
            if descendent not in subgraph.node_attr:
                subgraph.add_node(descendent)
                # print("adding subgraph node", descendent)
    for fname in subgraph.node_attr:
        for neighb in g.node_neighbors[fname]:
            if neighb in subgraph.node_attr:
                subgraph.add_edge(fname, neighb)
                # print("add subgraph edge", fname, neighb)
    return subgraph

def trim_inclusions_from_subgraph_extreme(
        fname,
        id,
        starting_nodes,
        orig_tg,
        current_tg,
        total_order,
        file_contents,
        dri,        
):
    if fname in dri.surrogates:
        compile_function = wrap_compile_w_surrogates(dri.surrogates[fname], id)
    else:
        builds, gold_objdump = generate_objdump_for_file(
            fname, id
        )
        if not builds:
            print("ERROR: could not compile", fname)
            return

        compile_function = wrap_compile_extreme(gold_objdump, id)

    prohibit_removal = prohibit_removal_of_non_subgraph_headers(
        dri, orig_tg, current_tg, starting_nodes)
    trim_unnecessary_headers_from_file(
        fname,
        current_tg,
        total_order,
        file_contents,
        compile_function,
        prohibit_removal,
        add_indirect_headers=False,
    )
    write_file(fname, file_contents[fname])


def remove_includes_from_subgraph_in_parallel(
        fnames,
        fork_manager,
        starting_nodes,
        orig_tg,
        current_tg,
        total_order,
        file_contents,
        dri
):
    for i, fname in enumerate(fnames):
        pid = fork_manager.mfork()
        if pid == 0:
            trim_inclusions_from_subgraph_extreme(
                fname,
                fork_manager.myjob_index+1,
                starting_nodes,
                orig_tg,
                current_tg,
                total_order,
                file_contents,
                dri
            )
            sys.exit(0)
    fork_manager.wait_for_remaining_jobs()

            
def remove_transcludes_from_subgraph(starting_files, fork_manager):
    dri = DontRemoveInclude()
    compilable_files, all_includes, file_contents, g = (
        create_inclusion_graph_from_source("pickled_include_graph.bin")
    )
    orig_tg = transitive_closure(g)

    add_inclusions_for_subgraph(orig_tg, starting_files, file_contents)
    subgraph = subgraph_from_starting_nodes(g, orig_tg, starting_files)
    print("created subgraph with", len(subgraph.node_attr), "nodes")
    equiv_sets = inclusion_equivalence_sets_from_desired_subgraph(subgraph)

    modified = False
    for count_round, es in enumerate(equiv_sets):
        equiv_set_filtered = []
        for fname in es:
            for starting_file in starting_files:
                if ( fname == starting_file or
                     starting_file in orig_tg.node_neighbors[fname] ):
                    print("  Examining file", fname)
                    equiv_set_filtered.append(fname)
                    break

        print("starting round", count_round, "with", len(equiv_set_filtered), "files to process")
        if len(equiv_set_filtered) == 0:
            print("No files to process in this round!")
            continue

        if count_round > 0 and modified:
            compilable_files, all_includes, file_contents, g = (
                create_inclusion_graph_from_source()
            )
            print("creating transitive closure")
            current_tg = transitive_closure(g)
            print("number of edges in transitive closure graph:", len(current_tg.edge_properties))
        else:
            current_tg = orig_tg
        total_order = total_order_from_graph(g)

        modified = True
        remove_includes_from_subgraph_in_parallel(
            equiv_set_filtered,
            fork_manager,
            starting_files,
            orig_tg,
            current_tg,
            total_order,
            file_contents,
            dri
        )

