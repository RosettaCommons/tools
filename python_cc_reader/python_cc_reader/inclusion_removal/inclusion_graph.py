# utility functions for working with
# inclusion graphs

import re, sys
import os
import pickle
from ..external.pygraph import pygraph
from toolz import curry


# from pygraph.algorithms.sorting import topological_sorting
from .add_headers import (
    add_autoheaders_to_file,
    add_autoheaders_to_filelines,
    write_file,
)
from .remove_header import remove_header_from_filelines
from .remove_duplicate_headers import remove_duplicate_headers_from_filelines
from .add_namespaces import (
    remove_using_namespace_from_lines,
    add_using_namespaces_to_lines,
    cleanup_auto_namespace_block_lines,
)
from .test_compile import (
    test_compile_from_lines,
    test_compile_from_stdin,
    generate_objdump_for_file,
    generate_objdump,
    test_compile_for_file_extreme,
    test_compile_extreme,
    test_compile_w_surrogates,
    cxxtest_test_compile,
)
# from ..cpp_parser.code_reader import (
#     heavy_code_reader,
# )
from ..cpp_parser.code_utilities import (
    known_circular_dependencies,
    scan_compilable_files,
    expand_includes_for_file,
    find_includes_at_global_scope,
    load_source_tree,
)
from .inclusion_equivalence_sets import inclusion_equivalence_sets
from . import dont_remove_include

# this is an optional way to get an inclusion graph into memory.
# the alternative is to get all the #includes from all
# reachable files in the source tree by calling
# scan_compilable_files in code_utilities.py or
# scan_files_to_create_inclusion_gaph
def read_inclusion_graph(filename):
    # the structure of the inclusion graph is simple:
    # all the nodes are listed first, one node per line
    # then all the edges are listed second, one edge per line, "B A" if B #includes A
    splitter = re.compile("\S+")
    file = open(filename, "r")
    done_with_nodes = False
    gr = pygraph.digraph()

    ignore = ["utility/keys/KeyLookup.functors.hh"]

    for line in file.readlines():
        toks = splitter.findall(line)
        if len(toks) == 0:
            continue
        if done_with_nodes:
            if len(toks) != 2:
                print(
                    "Error: While reading graph file",
                    filename,
                    ", encountered edge line",
                )
                print("Error:", line)
                print("Error: with", len(toks), "tokens (expected exactly 2 tokens)")
                sys.exit(1)
            else:
                if toks[0] not in gr.nodes():
                    print(
                        "Error: While reading graph file",
                        filename,
                        ", encountered edge line",
                    )
                    print("Error:", line)
                    print("Error: where the first node is not present in the graph")
                    sys.exit(1)
                if toks[1] not in gr.nodes():
                    print(
                        "Error: While reading graph file",
                        filename,
                        ", encountered edge line",
                    )
                    print("Error:", line)
                    print("Error: where the second node is not present in the graph")
                    sys.exit(1)
                if toks[1] in ignore:
                    continue
                gr.add_edge(toks[0], toks[1])
        else:
            if len(toks) == 2:
                done_with_nodes = True
                gr.add_edge(toks[0], toks[1])
            elif len(toks) == 1:
                if toks[0] in gr.nodes():
                    print("Error: found duplicate node", toks[0])
                    sys.exit(1)
                gr.add_node(toks[0])
            else:
                print(
                    "Error: While reading graph file",
                    filename,
                    ", encountered edge line",
                )
                print("Error:", line, end=" ")
                print(
                    "Error: with",
                    len(toks),
                    "tokens (expected either 1 token or 2 tokens)",
                )
                sys.exit(1)
    file.close()

    for node in gr.nodes():
        for neighb in gr.node_neighbors[node]:
            if node == neighb:
                print("Node connects to itself?!:", node)
                sys.exit(1)

    return gr


def write_inclusion_graph(g, filename):
    lines = []
    for node in g.nodes():
        lines.append(node + "\n")
    for node in g.nodes():
        for neighb in g.node_neighbors[node]:
            lines.append(node + " " + neighb + "\n")
    write_file(filename, lines)


def transitive_closure(g):
    retg = pygraph.digraph()
    retg.add_nodes(g.nodes())
    for node in retg.nodes():
        for neighb in g.node_neighbors[node]:
            retg.add_edge(node, neighb, label="direct")

    # input_dependents = {}
    # for node in input_dependents :
    #   input_dependents[ node ] = {}
    #   for neighb in g.node_neighbors[ node ] :
    #      input_dependents[ node ][ neighb ] = 0

    rtps = topological_sorting(retg)
    rtps.reverse()

    orig_neighb = {}
    reachable = {}
    for n in retg:
        reachable[n] = set(retg.node_neighbors[n])
        orig_neighb[n] = set(retg.node_neighbors[n])

    for n in rtps:
        for neighb in retg.node_neighbors[n]:
            reachable[n] |= reachable[neighb]
        for r in reachable[n]:
            if not r in orig_neighb[n]:
                retg.add_edge(n, r, label="indirect")

    return retg


def add_necessary_headers(tg, C_cc, headers):
    # Adds the headers to C_cc that are not included through transitive closure in C_cc already
    # and then returns the list of the headers that were added in topological sorted order
    topsort = topological_sorting(tg)
    added_headers = []
    for h in topsort:
        if h not in headers:
            continue
        if h in tg.node_neighbors[C_cc]:
            continue
        added_headers.append(h)
    add_autoheaders_to_file(C_cc, added_headers)
    return added_headers


# returns True when X_hh should not be removed from C_cc
def prohibit_removal_standard(C_cc, X_hh, DRI=None):
    is_ccfile = False
    is_header = False
    is_tmpl_header = False
    is_fwd_header = False
    is_hpp_header = False

    # don't touch anything in core/options
    if C_cc.find("core/options") != -1:
        return True

    if len(C_cc) < 3:
        return False

    term3 = C_cc[len(C_cc) - 3 :]
    if term3 == ".hh":
        is_header = True
    if term3 == ".cc":
        is_ccfile = True

    if len(C_cc) > 4:
        if C_cc[len(C_cc) - 4 :] == ".hpp":
            is_hpp_header = True

    if is_header:
        if len(C_cc) > 6:
            term7 = C_cc[len(C_cc) - 7 :]
            if term7 == ".fwd.hh":
                is_fwd_header = True
        if len(C_cc) > 7:
            term8 = C_cc[len(C_cc) - 8 :]
            if term8 == ".tmpl.hh":
                is_tmpl_header = True

    if is_fwd_header:
        if X_hh == "utility/pointer/owning_ptr.hh":
            return True
        if X_hh == "utility/pointer/access_ptr.hh":
            return True

    elif is_tmpl_header:
        source_header = C_cc[: len(C_cc) - 8] + ".hh"
        if X_hh == source_header:
            return True

    elif is_header:
        source_fwd_header = C_cc[: len(C_cc) - 3] + ".fwd.hh"
        # print source_fwd_header, X_hh
        if X_hh == source_fwd_header:
            return True
    elif is_hpp_header:
        source_fwd_header = C_cc[: len(C_cc) - 4] + ".fwd.hpp"
        if X_hh == source_fwd_header:
            return True
    elif is_ccfile:
        source_header = C_cc[: len(C_cc) - 3] + ".hh"
        if X_hh == source_header:
            return True
    else:
        # whoa -- what are we looking at?!
        print("Encountered weird file: ", C_cc)
        return True

    if DRI and DRI.leave_include_in_file(X_hh, C_cc):
        return True

    return False


def add_indirect_headers(tg, files):
    starts = [
        "utility/",
        "ObjexxFCL/",
        "numeric/",
        "basic/options/",
    ]  # skip any file that starts with this
    re_starts = {}
    for start in starts:
        re_starts[start] = re.compile(start)

    for fname in files:
        skip = False
        for start in starts:
            if re_starts[start].match(fname):
                skip = True
                break
        if skip:
            continue

        indirect_headers = []
        for neighb in tg.node_neighbors[fname]:
            if tg.edge_label(fname, neighb) == "indirect":
                indirect_headers.append(neighb)
        add_autoheaders_to_file(fname, indirect_headers)


# For files that are not part of the transitive-closure graph,
# detect the headers that have to be included and add those headers
# useful for getting the cxxtest headers straightened
def add_indirect_headers_to_file_outside_tg(tg, fname):
    direct_headers = find_includes_at_global_scope(fname, open(fname).readlines())
    direct_headers = set(direct_headers)
    indirect_headers = set()
    for header in direct_headers:
        if header not in tg.node_neighbors:
            continue  # not all inclusions are to files in tg
        for neighb in tg.node_neighbors[header]:
            if neighb not in indirect_headers and neighb not in direct_headers:
                indirect_headers.add(neighb)
    add_autoheaders_to_file(fname, indirect_headers)


# Remove any headers from the input file C_cc that are unnecessary for
# the compilation of C_cc.  C_cc may be a .cc, an .hh or a .fwd.hh file.
# compile_command should be a function that takes three arguments:
# C_cc, file_contents, and compile_args. It should return true if the
# compilation succeeds and false if the compilation fails
def trim_unnecessary_headers_from_file(
        C_cc,
        tg,
        total_order,
        file_contents,
        compile_command, 
        prohibit_removal,
        add_indirect_headers=True,
        super_cautious=False,
):
    """
    compiled_command: (possibly curried) function: arg1: fname,
                      arg2: dict w/ file lines
    prohibit_removal: (possibly curried) function: arg1: fname,
                      arg2: included header w/i arg1 file
    """
    print("Examining includes within", C_cc)
    sys.stdout.flush()

    if not compile_command(C_cc, file_contents):
        print("Skipping ", C_cc, " since it does not compile as is")
        return

    backup_C_cc = file_contents[C_cc][:]  # deep copy
    last_compiling_version = backup_C_cc[:]

    file_contents[C_cc] = remove_duplicate_headers_from_filelines(
        C_cc, file_contents[C_cc]
    )

    C_cc_headers = None
    indirect_headers = []
    
    if C_cc in tg.node_neighbors:
        # if we're dealing with a file that's part of the tg,
        # then we can use the edges in tg to construct the set
        # of indirect_headers
        C_cc_headers = tg.node_neighbors[C_cc]
        for header in C_cc_headers:
            # print header, tg.edge_label( C_cc, header ), C_cc
            if tg.edge_label(C_cc, header) == "indirect":
                indirect_headers.append(header)
    else:
        # if we're not dealing with a file that's part of the
        # tg, then read the headers from the file_contents
        # add in transitive headers for any files that are
        # part of the tg.  (Note there may be headers which
        # are not part of tg). Note that if we're going to
        # remove includes for a file, then that file must be
        # included in the file_contents map, even though it
        # doesn't have to be part of the tg
        assert C_cc in file_contents
        C_cc_headers = find_includes_at_global_scope(C_cc, file_contents[C_cc])
        if len(C_cc_headers) == 0:
            print("FOUND NO HEADERS AT GLOBAL SCOPE!")
        C_cc_headers_set = set(C_cc_headers)
        indirect_set = set()
        for header in C_cc_headers:
            if header in tg:
                for neighb in tg.node_neighbors[header]:
                    indirect_set.add(neighb)
        indirect_headers = list(indirect_set)

    for h in C_cc_headers:  # temp debug
        print("File", C_cc, " has header:", h)  # temp debug
    for ih in indirect_headers:  # temp debug
        print("File", C_cc, " has indirect header:", ih)  # temp debug

    if add_indirect_headers:
        print("Adding auto-header to", C_cc)
        file_contents[C_cc] = add_autoheaders_to_filelines(
            C_cc, file_contents[C_cc], indirect_headers
        )

        if compile_command(C_cc, file_contents):
            print("autoheader addition succeeded")
        else:
            print("Warning: could not resolve compilation issues with ", C_cc)
            file_contents[C_cc] = backup_C_cc
            return

    backup_C_cc = file_contents[C_cc][:]  # deep copy
    last_compiling_version = backup_C_cc[:]

    print("Adding using namespaces to lines")
    namespaces, file_contents[C_cc] = add_using_namespaces_to_lines(
        C_cc, file_contents[C_cc]
    )

    if not compile_command(C_cc, file_contents):
        print(
            "Warning: failed to compile",
            C_cc,
            "after adding using namespace declarations",
        )
        file_contents[C_cc] = backup_C_cc
        newnamespaces = []
        for ns in namespaces:
            print("Testing the addition of namespace", ns)
            backup_C_cc = file_contents[C_cc][:]  # deep copy
            test_namespaces = list(newnamespaces)
            test_namespaces.append(ns)
            blah, file_contents[C_cc] = add_using_namespaces_to_lines(
                C_cc, file_contents[C_cc], test_namespaces
            )
            if not compile_command(C_cc, file_contents):
                file_contents[C_cc] = backup_C_cc
            else:
                last_compiling_version = file_contents[C_cc][:]
                newnamespaces.append(ns)
        namespaces = newnamespaces

    for namespace in namespaces:
        print(namespace, " added in auto-namespace block")

    if not compile_command(C_cc, file_contents):
        print("ERROR: failed to compile after namespace addition twice")
        file_contents[C_cc] = last_compiling_version
        return

    C_cc_headers_in_tg = []
    for X_hh in C_cc_headers:
        if X_hh in tg:
            C_cc_headers_in_tg.append(X_hh)
        else:
            print("Header", X_hh, "not in tg. Skipping.")

    headers = topologically_sorted(total_order, C_cc_headers_in_tg)
    necessary = []  # tg.node_incidence[ C_cc ]
    for X_hh in headers:

        print("...Testing: is ", X_hh, "necessary for", C_cc, end=" ")

        if prohibit_removal(C_cc, X_hh):
            print("yes.  Prohibited from removing", X_hh)
            necessary.append(X_hh)
            continue

        X_hh_necessary = False
        for nec in necessary:
            if X_hh in tg.node_neighbors[nec]:
                X_hh_necessary = True
                if X_hh in tg:
                    print(
                        "yes.  Dependent of",
                        X_hh,
                        "(",
                        len(tg.node_neighbors[X_hh]),
                        ") is necessary.",
                    )
                else:
                    print("yes.  Dependent of", X_hh, "is necessary.")
                break

        if X_hh_necessary:
            # only remove the header if it's not an auto-header
            # by adding the fourth argument (True) here, this directs
            # the remove_header_from_filelines function to leave
            # non-auto-headers intact.
            if super_cautious:
                contents_after_removal = remove_header_from_filelines(
                    C_cc, file_contents[C_cc], X_hh, True
                )

                # OK super extra slow.  Make sure that removing this header didn't break the build.
                if file_contents[C_cc] != contents_after_removal:
                    backup_C_cc = file_contents[C_cc][:]  # deep copy
                    file_contents[C_cc] = contents_after_removal
                    if not compile_command(C_cc, file_contents):
                        print(
                            "Whoa!  Removing header ",
                            X_hh,
                            " which is transitively included breaks the build! WTF?!",
                        )
                        file_contents[C_cc] = backup_C_cc
                    else:
                        print(
                            "Safely removed ",
                            X_hh,
                            "which is transitively included and an auto-header",
                        )
            else:
                file_contents[C_cc] = remove_header_from_filelines(
                    C_cc, file_contents[C_cc], X_hh, True
                )

            continue

        backup_C_cc = file_contents[C_cc][:]  # deep copy
        file_contents[C_cc] = remove_header_from_filelines(
            C_cc, file_contents[C_cc], X_hh
        )
        if not compile_command(C_cc, file_contents):
            if X_hh in tg:
                print(
                    "yes.  Removing",
                    X_hh,
                    "(",
                    len(tg.node_neighbors[X_hh]),
                    ") breaks the build.",
                )
            else:
                print("yes.  Removing", X_hh, "breaks the build.")
            file_contents[C_cc] = backup_C_cc
            necessary.append(X_hh)
        else:
            last_compiling_version = file_contents[C_cc][:]
            if X_hh in tg:
                print(
                    "no.",
                    X_hh,
                    "(",
                    len(tg.node_neighbors[X_hh]),
                    ") may be safely removed.",
                )
            else:
                print("no.", X_hh, "may be safely removed.")

    for ns in namespaces:
        backup_C_cc = file_contents[C_cc][:]  # deep copy
        file_contents[C_cc] = remove_using_namespace_from_lines(
            C_cc, file_contents[C_cc], ns
        )
        print("...Testing: is namespace", ns, "necessary?", end=" ")
        if not compile_command(C_cc, file_contents):
            # if not test_compile_from_lines( expand_includes_for_file( C_cc, file_contents ) ) :
            print("yes. using namespace", ns, "must remain")
            file_contents[C_cc] = backup_C_cc
        else:
            last_compiling_version = file_contents[C_cc][:]
            print("no.  Namespace", ns, " does not need a using declaration")

    # remove the auto-namespace block if it's empty
    file_contents[C_cc] = cleanup_auto_namespace_block_lines(file_contents[C_cc])

    if not compile_command(C_cc, file_contents):
        print(
            "Warning: Could not compile",
            C_cc,
            "at conclusion of trim_unnecessary_headers_from_file",
        )
        open("blah", "w").writelines(file_contents[C_cc])
        file_contents[C_cc] = last_compiling_version
        if not super_cautious:
            # try again with super-slow removal
            trim_unnecessary_headers_from_file(
                C_cc,
                tg,
                total_order,
                file_contents,
                compile_command,
                prohibit_removal,
                add_indirect_headers,
                super_cautious=True,
            )

    print("Finished Examining includes within", C_cc)
    sys.stdout.flush()

# this function takes a list of files to process.  It then reads the entire
# source tree, computes the transitive closure graph, a total-order for
# the graph, and then finally begins to trim headers from the individual files
# that it has been asked to process.  After processing each file, it then
# writes that file to disk.
def trim_inclusions_from_files(filelist):
    print("loading source tree")
    compilable_files, all_includes, file_contents = load_source_tree()
    print("creating graph from includes")
    g = create_graph_from_includes(all_includes)
    remove_known_circular_dependencies_from_graph(g)
    print("creating transitive closure")
    tg = transitive_closure(g)

    print("computing total order")
    total_order = total_order_from_graph(g)
    dummy = 1
    dri = dont_remove_include.DontRemoveInclude()
    for fname in filelist:
        print("beginning work on", fname)
        trim_unnecessary_headers_from_file(
            fname, tg, total_order, file_contents, simple_compile_from_lines,
            canonical_prohibit_removal(dri)
        )
        write_file(fname, file_contents[fname])


# Similar to trim_inclusions_from_files, except it relies on the
# "extreme" version of the compilation test which also compares
# the objdump output from the compiled object file
def trim_inclusions_from_files_extreme(
        filelist,
        id,
        file_contents,
        tg,
        total_order,
        dri,
        super_cautious=False):

    for fname in filelist:
        trim_inclusions_from_file_extreme(
            fname,
            id,
            file_contents,
            tg,
            total_order,
            dri,
            super_cautious)

        
def trim_inclusions_from_file_extreme(
        fname,
        id,
        file_contents,
        tg,
        total_order,
        dri,
        super_cautious=False):

    if fname in dri.surrogates:
        compile_function = wrap_compile_w_surrogates(dri.surrogates[fname], id)
    else:
        builds, gold_objdump = generate_objdump_for_file(
            fname, id
        )
        if not builds:
            print("ERROR: could not compile", fname)
            return
        else:
            compile_function = wrap_compile_extreme(gold_objdump, id)

    trim_unnecessary_headers_from_file(
        fname,
        tg,
        total_order,
        file_contents,
        compile_function,
        canonical_prohibit_removal(dri),
        add_indirect_headers=False,
        super_cautious=super_cautious,
    )
    write_file(fname, file_contents[fname])



# Similar to trim_inclusions_from_files, except it relies on the
# "extreme" version of the compilation test which also compares
# the objdump output from the compiled object file
def trim_inclusions_from_cxxtest_files(filelist, id):
    compilable_files, all_includes, file_contents = load_source_tree()
    g = create_graph_from_includes(all_includes)
    remove_known_circular_dependencies_from_graph(g)
    tg = transitive_closure(g)
    total_order = total_order_from_graph(g)
    dri = dont_remove_include.DontRemoveInclude()
    for fname in filelist:
        file_contents[fname] = open(fname).readlines()
        trim_unnecessary_headers_from_file(
            fname, tg, total_order, file_contents, wrap_cxxtest_compile(id),
            canonical_prohibit_removal(dri)
        )
        write_file(fname, file_contents[fname])

@curry
def wrap_compile_extreme(gold_objdump, id, fname, file_contents):
    write_file(fname, file_contents[fname])
    return test_compile_for_file_extreme(
        fname, gold_objdump, id
    )

@curry
def wrap_compile_w_surrogates(surrogates, id, fname, file_contents):
    write_file(fname, file_contents[fname])
    return test_compile_w_surrogates(fname, surrogates, id)


def simple_compile_from_lines(C_cc, file_contents):
    return test_compile_from_stdin(C_cc, file_contents)

@curry
def wrap_cxxtest_compile(id, C_cc, file_contents):
    # no option to read from stdin, so first write the file to disk
    open(C_cc, "w").writelines(file_contents[C_cc])
    return cxxtest_test_compile(C_cc, False, id)


def backup_file(filename, extension):
    newlines = open(filename).readlines()
    write_file(filename + "." + extension, newlines)


def restore_backup(filename, extension):
    newlines = open(filename + "." + extension).readlines()
    write_file(filename, newlines)


def create_graph_from_includes(all_includes):
    g = pygraph.digraph()
    for fname in list(all_includes.keys()):
        if fname not in g.nodes():
            g.add_node(fname)
        for include in all_includes[fname]:
            if include not in g.nodes():
                g.add_node(include)
            g.add_edge(fname, include)
    return g


def scan_files_to_create_inclusion_graph():
    all_includes = scan_compilable_files()
    return create_graph_from_includes(all_includes)


def remove_known_circular_dependencies_from_graph(dg):
    circular_dependencies = known_circular_dependencies()
    for cdep in circular_dependencies:
        if dg.has_edge(cdep[0], cdep[1]):
            dg.del_edge(cdep[0], cdep[1])


def total_order_from_graph(dg):
    toposort = topological_sorting(dg)
    total_order = {}
    count_total_order = 0
    for file in toposort:
        count_total_order += 1
        total_order[file] = count_total_order
    return total_order


# Return the contents of file_list in topological-sorted-order, given
# a total order
def topologically_sorted(total_order, file_list):
    headers = []
    tot_order_set = []
    for fname in file_list:
        tot_order_set.append((fname, total_order[fname]))
    values = sorted(tot_order_set, key = lambda x: x[1])
    # print "Dependent headers for ", C_cc
    for val in values:
        # print "   ", val[ 0 ], val[ 1 ]
        headers.append(val[0])
    return headers


# This code is a simple extension of the find_cycle subroutine of
# the pygraph library.
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

    if type(graph) == pygraph.digraph:
        directed = True
    else:
        raise InvalidGraphType

    def find_cycle_to_ancestor(node, ancestor):
        """
        Find a cycle containing both node and ancestor.
        """
        path = []
        orignode = node
        while node != ancestor:
            if node is None:
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
            if cycle:
                return
            if each not in visited:
                spanning_tree[each] = node
                dfs(each)
            else:
                if directed or spanning_tree[node] is not each:
                    cycle.extend(find_cycle_to_ancestor(node, each))

    visited = {}  # List for marking visited and non-visited nodes
    spanning_tree = {}  # Spanning tree
    cycle = []
    cycles = []

    # Algorithm outer-loop
    for each in graph:
        # Select a non-visited node
        if each not in visited:
            spanning_tree[each] = None
            root = each
            # Explore node's connected component
            dfs(each)
            if cycle:
                cycles.append(cycle)
                cycle = []

    return cycles


def print_cycles_in_source_tree():

    gr = scan_files_to_create_inclusion_graph()
    remove_known_circular_dependencies_from_graph(gr)

    cycles = find_cycles(gr)

    if len(cycles) == 0:
        print("Found no cycles in inclusion graph")
    else:

        count = 0
        for cycle in cycles:
            count += 1
            print("-------")
            print("Cycle #", count)
            for node in cycle:
                print(node)


def desired_node(node):
    name_list = node.split("/")
    if len(name_list) == 1:
        return False
    desired_libs = ["utility", "numeric", "basic", "core", "protocols", "devel", "apps"]
    if name_list[0] in desired_libs:
        return True
    return False


def inclusion_equivalence_sets_from_desired_subgraph(g):
    gprime = pygraph.digraph()
    for node in g.nodes():
        if desired_node(node):
            gprime.add_node(node)
    for edge in g.edges():
        if desired_node(edge[0]) and desired_node(edge[1]):
            gprime.add_edge(edge[0], edge[1])
    equiv_sets = inclusion_equivalence_sets(gprime)
    return equiv_sets


def write_equiv_sets_file(fname, equiv_sets):
    lines = []
    for es in equiv_sets:
        for include in es:
            lines.append(include + "\n")
        lines.append("\n")  # blank line to separate equiv_sets
    write_file(fname, lines)


def read_equiv_sets_file(fname):
    equiv_sets = []
    lines = open(fname).readlines()
    es = []
    for line in lines:
        if line == "\n":
            equiv_sets.append(es)
            es = []
        else:
            es.append(line.strip())
    return equiv_sets


#### CODE STOLEN FROM PYGRAPH LIBRARY
# Copyright (c) 2007-2009 Pedro Matiello <pmatiello@gmail.com>
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation
# files (the "Software"), to deal in the Software without
# restriction, including without limitation the rights to use,
# copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following
# conditions:

# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.


"""
Search algorithms.

@sort: breadth_first_search, depth_first_search
"""


# Depth-first search


def depth_first_search(graph):
    """
    Depth-first search.

    @type  graph: graph, digraph
    @param graph: Graph.
    
    @type  root: node
    @param root: Optional root node (will explore only root's connected component)

    @rtype:  tuple
    @return: A tupple containing a dictionary and two lists:
        1. Generated spanning tree
        2. Graph's preordering
        3. Graph's postordering
    """

    def dfs(node):
        """
        Depht-first search subfunction.
        """
        visited[node] = 1
        pre.append(node)
        # Explore recursively the connected component
        for each in graph[node]:
            if each not in visited:
                spanning_tree[each] = node
                dfs(each)
        post.append(node)

    visited = {}  # List for marking visited and non-visited nodes
    spanning_tree = {}  # Spanning tree
    pre = []  # Graph's preordering
    post = []  # Graph's postordering

    # Algorithm loop
    for each in graph:
        # Select a non-visited node
        if each not in visited:
            spanning_tree[each] = None
            # Explore node's connected component
            dfs(each)

    return spanning_tree, pre, post


# Topological sorting
def topological_sorting(graph):
    """
    Topological sorting.

    @attention: Topological sorting is meaningful only for directed acyclic graphs.

    @type  graph: digraph
    @param graph: Graph.

    @rtype:  list
    @return: Topological sorting for the graph.
    """
    # The topological sorting of a DAG is equivalent to its reverse postordering.
    tmp, tmp, post = depth_first_search(graph)
    post.reverse()
    return post

@curry
def canonical_prohibit_removal(DRI, C_cc, X_hh):
    """ Only prohibit the removal of the standard
    set of headers from C_cc
    """
    return prohibit_removal_standard(C_cc, X_hh, DRI)

@curry
def prohibit_removal_of_non_subgraph_headers(
        DRI,
        original_tg,
        current_tg,
        starting_nodes,
        C_cc,
        X_hh
):
    """Prohibit the removal of headers that are not transcluded by
    a file which C_cc itself transcludes"""
    if prohibit_removal_standard(C_cc, X_hh, DRI): 
       return True
    if X_hh in starting_nodes:
        return False
    if C_cc in starting_nodes:
        return False

    for starting_node in starting_nodes:
        # print("starting node", starting_node)
        # print("C_cc", C_cc)
        # print("neighbors", current_tg.node_neighbors[C_cc])
        # print("starting_node in current_tg.node_neighbors[C_cc]", starting_node in current_tg.node_neighbors[C_cc])
        # print("X_hh in original_tg.node_neighbors[starting_node]", X_hh in original_tg.node_neighbors[starting_node])
        if (
            starting_node in current_tg.node_neighbors[C_cc] and
            X_hh in original_tg.node_neighbors[starting_node]
        ):
            return False
    print("file", X_hh, "is not a trasclude of any of the starting nodes; leave it...", end="")
    # Debugging code
    # if X_hh == "utility/CSI_Sequence.hh":
    #     for starting_node in starting_nodes:
    #         print("starting node", starting_node)
    #         print("C_cc", C_cc)
    #         print("C_cc neighbors", current_tg.node_neighbors[C_cc])
    #         print("starting node neighbors", original_tg.node_neighbors[starting_node])
    #         print("starting_node in current_tg.node_neighbors[C_cc]", starting_node in current_tg.node_neighbors[C_cc])
    #         print("X_hh in original_tg.node_neighbors[starting_node]", X_hh in original_tg.node_neighbors[starting_node])
    #     print("FAILED TO DETECT REMOVABLE HEADER!!!")
    #     sys.exit(0)
    return True

# report the number of #includes for each file.
if __name__ == "__main__":
    g = scan_files_to_create_inclusion_graph()
    for node in g.nodes():
        print(node, len(g.node_neighbors[node]))
