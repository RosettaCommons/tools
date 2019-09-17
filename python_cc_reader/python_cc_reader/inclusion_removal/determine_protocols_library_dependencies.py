import sys

from ..cpp_parser import code_utilities
from . import library_levels
from ..external.pygraph import pygraph
from . import inclusion_graph


# $1 == --graphviz
# $2 == target prefix name

# input:  the inclusion graph of the entire source tree
# output: the graph describing the dependencies between the various protocols libraries
def create_protocols_interlibrary_dependency_graph(g):
    prot_levels = library_levels.protocols_levels()

    protocols_library_graph = pygraph.digraph()
    projects = code_utilities.rosetta_projects()
    for lib in projects["src"]:
        # print lib
        if lib.find("protocols") != -1:
            suffix = lib[9:]
            if len(suffix) == 2:
                protocols_library_graph.add_node(lib[-1:])  # e.g. "7" for protocols.7
            elif len(suffix) == 4:
                protocols_library_graph.add_node(
                    lib[-1:] + lib[-3:-2]
                )  # e.g. 4a for protocols_a.4

    for edge in g.edges():
        file1, file2 = edge
        lib1 = library_levels.libname_for_file(file1)
        lib2 = library_levels.libname_for_file(file2)
        if lib1 != "protocols" or lib2 != "protocols":
            continue
        lev1, col1 = library_levels.library_and_column_for_file(prot_levels, file1)
        lev2, col2 = library_levels.library_and_column_for_file(prot_levels, file2)
        node1 = str(lev1) + col1
        node2 = str(lev2) + col2
        if node2 not in protocols_library_graph.neighbors(node1):
            # print "adding edge", node1, node2
            protocols_library_graph.add_edge(node1, node2)

    return protocols_library_graph


# input:  the inclusion graph of the entire source tree
# output: a set of lines for a dot-graph describing the interlibrary dependencies within prorocols
def create_protocols_interlibrary_dependency_dot_file(g):
    prot_levels = library_levels.protocols_levels()
    protocols_library_graph = create_protocols_interlibrary_dependency_graph(g)
    outlines = []
    endl = "\n"
    outlines.append("digraph G {" + endl)
    for level in range(len(prot_levels)):
        line = "subgraph level" + str(level + 1) + " { rank = same; "
        if len(prot_levels[level]) == 1:
            line += '"' + str(level + 1) + '";}' + endl
        else:
            for column in range(len(prot_levels[level])):
                line += '"' + str(level + 1) + chr(ord("a") + column) + '";'
            line += "}" + endl
        outlines.append(line)
    for node in protocols_library_graph.nodes():
        outlines.append('"' + node + '";' + endl)
    for node in protocols_library_graph.nodes():
        for node2 in protocols_library_graph.neighbors(node):
            if node == node2:
                continue  # skip loops
            outlines.append('"' + node2 + '" -> "' + node + '" [ dir=back ];' + endl)
    outlines.append("}" + endl)
    return outlines


# input: the name of a particular .src.settings file
# output: the library prefix and a set of all the toplevel directories within that file
def toplevel_directories_from_srcsettings_file(srcsettings_fname):
    f = open(srcsettings_fname).read()
    exec(f)
    topleveldirs = set([])
    libname = ""
    for dirname in sources:
        cols = dirname.split("/")
        if libname == "":
            libname = cols[0]
        else:
            assert libname == cols[0]
        if len(cols) < 2:
            continue
        if cols[1] not in topleveldirs:
            topleveldirs.add(cols[1])
    return libname, topleveldirs


# input:  the inclusion graph for the entire source tree, a particular .src.settings file name defining one library
# output: a graph describing the relationship between nodes in a particular library
def create_intralibrary_dependency_graph(g, srcsettings_fname):
    libname, topleveldirs = toplevel_directories_from_srcsettings_file(
        srcsettings_fname
    )
    topleveldir_graph = pygraph.digraph()
    for topleveldir in topleveldirs:
        topleveldir_graph.add_node(topleveldir)

    liblevel, libcolumn = library_level_and_column_from_src_settingsfile_name(
        libname, srcsettings_fname
    )
    prot_levels = library_levels.protocols_levels()

    for edge in g.edges():
        file1, file2 = edge
        lib1 = library_levels.libname_for_file(file1)
        lib2 = library_levels.libname_for_file(file2)
        if lib1 != libname or lib2 != libname:
            continue

        lev1, col1 = library_levels.library_and_column_for_file(prot_levels, file1)
        if lev1 != liblevel or col1 != libcolumn:
            continue
        lev2, col2 = library_levels.library_and_column_for_file(prot_levels, file2)
        if lev2 != liblevel or col2 != libcolumn:
            continue
        node1 = file1.split("/")[1]
        node2 = file2.split("/")[1]
        if node1 not in topleveldirs:
            print(
                "ERROR: node1:",
                node1,
                "not in the topleveldirs set in spite of being identified by library_and_column_for_file as belonging to the target library",
            )
            print("all nodes:", topleveldirs)
            assert node1 in topleveldirs
        if node2 not in topleveldirs:
            print(
                "ERROR: node2:",
                node2,
                "not in the topleveldirs set in spite of being identified by library_and_column_for_file as belonging to the target library",
            )
            print("all nodes:", topleveldirs)
            assert node2 in topleveldirs
        if node2 not in topleveldir_graph.neighbors(node1):
            # print "adding edge", node1, node2
            topleveldir_graph.add_edge(node1, node2)

    return topleveldir_graph


def create_intralibrary_dependency_dot_file(g, srcsettings_fname):
    topleveldir_graph = create_intralibrary_dependency_graph(g, srcsettings_fname)
    lines = []
    lines.append("digraph G { rankdir=BT;\n")
    for node in topleveldir_graph.nodes():
        lines.append('"' + node + '";\n')
        for neighbor in topleveldir_graph.neighbors(node):
            if node == neighbor:
                continue
            lines.append('"' + node + '" -> "' + neighbor + '";\n')
    lines.append("}\n")
    return lines


# returns a tuple -- an integer and a string -- representing the library name and column.
def library_level_and_column_from_src_settingsfile_name(
    libprefix, lib_srcsettings_fname
):
    assert lib_srcsettings_fname.find(libprefix) != -1
    suffix = lib_srcsettings_fname[len(libprefix) :][:-13]
    # print libprefix, suffix, lib_srcsettings_fname[len(libprefix):]
    if len(suffix) == 2:
        return int(suffix[-1:]), ""
    elif len(suffix) == 4:
        return int(suffix[-1:]), suffix[-3:2]


def draw_all_protocol_dependencies(g):
    colors_for_level = [
        "black",
        "indigo",
        "crimson",
        "green",
        "chocolate",
        "cyan3",
        "blue1",
        "salmon2",
    ]

    intralibrary_dependency_graphs = {}
    projects = code_utilities.rosetta_projects()
    for lib in projects["src"]:
        if lib.find("protocols") == -1:
            continue
        libsrcsettings = lib + ".src.settings"
        level, column = library_level_and_column_from_src_settingsfile_name(
            "protocols", libsrcsettings
        )
        libid = str(level) + column
        intralibrary_dependency_graphs[libid] = create_intralibrary_dependency_graph(
            g, libsrcsettings
        )
    prot_levels = library_levels.protocols_levels()
    protocols_library_graph = create_protocols_interlibrary_dependency_graph(g)

    endl = "\n"
    graph_lines = []
    graph_lines.append("digraph G { compound=true; rankdir=BT; " + endl)
    representative_node = {}

    for level in range(len(prot_levels)):
        # ranksameline = "{ rank=same; "
        # graph_lines.append( "subgraph cluster_level" + str(level+1) + "{\n" )
        for column in range(len(prot_levels[level])):
            if len(prot_levels[level]) == 1:
                sublibname = str(level + 1)
            else:
                sublibname = str(level + 1) + chr(ord("a") + column)
            line = (
                'subgraph "cluster_'
                + sublibname
                + '" { label="protocols.'
                + sublibname
                + '"; color='
                + colors_for_level[int(sublibname[:1]) - 1]
                + ";"
            )
            intralib_dep_graph = intralibrary_dependency_graphs[sublibname]
            representative_node[sublibname] = intralib_dep_graph.nodes()[0]
            # ranksameline += intralib_dep_graph.nodes()[ 0 ] + " "
            for node in intralib_dep_graph.nodes():
                line += ' "' + node + '"; '
                for neighbor in intralib_dep_graph.neighbors(node):
                    if node == neighbor:
                        continue  # skip loops
                    line += ' "' + node + '" -> "' + neighbor + '";'
            line += "}" + endl
            graph_lines.append(line)
        # if len(prot_levels[level]) > 1 :
        #   for column in xrange(1, len(prot_levels[level])) :
        #      sublibname = str(level+1) + chr( ord('a') + column )
        #      prevsublibname = str(level+1) + chr( ord('a') + column-1 )
        #      graph_lines.append( '{rank=same; "' + representative_node[ sublibname ] + '"; "' + representative_node[ prevsublibname ] + '"}\n' )
        # graph_lines.append( "}\n" )

        # ranksameline += "}\n"
        # graph_lines.append( ranksameline )

    # now draw dependencies between the libraries
    for node in protocols_library_graph.nodes():
        nodelevel = int(node[:1])
        nodecolor = colors_for_level[nodelevel - 1]
        for node2 in protocols_library_graph.neighbors(node):
            if node == node2:
                continue  # skip loops
            graph_lines.append(
                '"'
                + representative_node[node]
                + '" -> "'
                + representative_node[node2]
                + '" [ minlen=10, lhead="cluster_'
                + node2
                + '", ltail="cluster_'
                + node
                + '", color='
                + nodecolor
                + " ];"
                + endl
            )
            # for n1 in intralibrary_dependency_graphs[ node ].nodes() :
            #   for n2 in intralibrary_dependency_graphs[ node2 ].nodes() :
            #      graph_lines.append( '"' + n1 + '" -> "' + n2 + '" [ style=invis ] \n' )
    graph_lines.append("}" + endl)
    return graph_lines


if __name__ == "__main__":

    graphviz_output = False
    if len(sys.argv) > 2 and sys.argv[1] == "--graphviz":
        graphviz_output = True
        output_prefix = sys.argv[2]

    compilable_files, all_includes, file_contents = code_utilities.load_source_tree()
    print("...source tree loaded")
    g = inclusion_graph.create_graph_from_includes(all_includes)
    print("...inclusion graph created")

    # protocols_tg = inclusion_graph.transitive_closure( protocols_library_graph )
    if not graphviz_output:
        protocols_tg = protocols_library_graph
        for node in protocols_tg.nodes():
            print("Dependencies for node", node)
            for node_neighbor in protocols_tg.neighbors(node):
                print("  ", node_neighbor)
    else:
        graphlines = draw_all_protocol_dependencies(g)
        open(output_prefix + "_all_dependencies.dot", "w").writelines(graphlines)

        graphlines = create_protocols_interlibrary_dependency_dot_file(g)
        open(output_prefix + "_interlibrary_dependencies.dot", "w").writelines(
            graphlines
        )

        projects = code_utilities.rosetta_projects()
        for lib in projects["src"]:
            if lib.find("protocols") == -1:
                continue
            libsrcsettings = lib + ".src.settings"
            graphlines = create_intralibrary_dependency_dot_file(g, libsrcsettings)
            open(
                output_prefix + "_" + lib + "_intralib_dependencies.dot", "w"
            ).writelines(graphlines)
