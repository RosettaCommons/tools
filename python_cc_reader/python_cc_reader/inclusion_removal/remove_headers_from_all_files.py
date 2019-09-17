import subprocess
import math
import random
import re
import sys

from ..cpp_parser.code_utilities import load_source_tree, expand_includes_for_file

from .inclusion_graph import create_graph_from_includes, remove_known_circular_dependencies_from_graph, transitive_closure, add_indirect_headers, inclusion_equivalence_sets_from_desired_subgraph, write_equiv_sets_file, trim_inclusions_from_files_extreme
from .test_compile import test_compile, tar_everything
#from .inclusion_equivalence_sets import *
#from .add_headers import *
#from .add_namespaces import *
#from .remove_header import 
#from .remove_duplicate_headers import *
#from ..cpp_parser.code_reader import *
#from .reinterpret_objdump import *

from .dont_remove_include import DontRemoveInclude
#from ...cpp_parser import code_reader
#from ...external.pygraph import pygraph

from ..utility.fork_manager import ForkManager


def whole_shebang(fork_manager: ForkManager, skip_test_compile=False):

    # OVERVIEW:
    # first create the include graph
    # trim this graph to remove known cycles (vector1_bool, vector0_bool)

    # second, verify that all files compile on their own

    # third, compute the transitive closure graph

    # fourth, add all headers in transitive closure

    # fifth, compute the equivalence sets, focusing only on files in
    # core/ protocols/ devel/ and apps/

    # sixth, iterate across all equivalence sets
    #   divide up the equivalence sets among the compilable files
    #   distribute tasks
    #   each node reads the source tree, and computes the transitive closure
    #   subgraph, looking only at the files proceeded by the files covered by
    #   this task* (*write this subroutine)
    #   Then
    #   iterating across files in this subset:
    #      if the file contains nothing but headers to other headers, continue
    #      remove all inessential headers from the file
    #      through successive recompiles
    #      write out the finished file

    # BEGIN CODE

    # first create the include graph
    # trim this graph to remove known cycles (vector1_bool, vector0_bool)

    compilable_files, all_includes, file_contents = load_source_tree()
    g = create_graph_from_includes(all_includes)
    remove_known_circular_dependencies_from_graph(g)

    # second, verify that all files compile on their own
    any_fail_to_compile = False
    def set_any_failed(_1, _2):
        any_fail_to_compile
    fork_manager.error_callback = set_any_failed
    
    if not skip_test_compile:
        for fname in compilable_files :
            print("testing compilation of", fname)
            pid = fork_manager.mfork()
            if pid == 0:
                if not test_compile(fname) :
                    print( "Error: ", fname, "does not compile on its own")
                    sys.exit(1)
                else:
                    sys.exit(0)
    fork_manager.wait_for_remaining_jobs()
    if any_fail_to_compile :
        print("Error: coud not compile all files on their own")
        sys.exit(0)

    # third, compute the transitive closure graph
    tg = transitive_closure(g)

    # fourth, add all headers in transitive closure
    tar_everything("bu_starting_code")
    add_indirect_headers(tg, compilable_files)
    tar_everything("bu_transclose_headers_round_0")

    # fifth, compute the equivalence sets, focusing only on files in
    # core/ protocols/ devel/ and apps/
    equiv_sets = inclusion_equivalence_sets_from_desired_subgraph(g)
    write_equiv_sets_file("equivalence_sets_wholeshebang.txt", equiv_sets)

    # equiv_sets = read_equiv_sets_file( "equivalence_sets_wholeshebang.txt" )

    #funcs = (
    #    test_compile_from_lines,
    #    central_compile_command,
    #    remove_duplicate_headers_from_filelines,
    #    no_empty_args,
    #    test_compile_extreme,
    #    generate_objdump,
    #    labeled_instructions,
    #    ignore_instructions,
    #    regexes_for_instructions,
    #    skip_sections,
    #    relabel_sections,
    #    compare_objdump_lines,
    #    add_autoheaders_to_filelines,
    #    remove_header_from_filelines,
    #    group_and_sort_headers,
    #    group_headers,
    #    add_using_namespaces_to_lines,
    #    remove_using_namespace_from_lines,
    #    cleanup_auto_namespace_block_lines,
    #    expand_includes_for_file,
    #    follow_includes_for_file,
    #    auto_ns_comment,
    #    transitive_closure,
    #    prohibit_removal,
    #    trim_unnecessary_headers_from_file,
    #    trim_inclusions_from_files_extreme,
    #    wrap_compile_extreme,
    #    load_source_tree,
    #    create_graph_from_includes,
    #    remove_known_circular_dependencies_from_graph,
    #    known_circular_dependencies,
    #    topological_sorting,
    #    depth_first_search,
    #    total_order_from_graph,
    #    transitive_closure,
    #    write_file,
    #    topologically_sorted,
    #    libraries_with_ccfiles_to_examine,
    #    scan_compilable_files,
    #    libraries_with_hhfiles_to_examine,
    #    directories_with_ccfiles_to_examine,
    #    directories_with_hhfiles_to_examine,
    #    include_for_line,
    #    find_all_includes,
    #    find_includes_at_global_scope,
    #    compiled_cc_files,
    #    strip_toendofline_comment,
    #)

    # modules = (
    #     "re",
    #     "subprocess",
    #     "code_reader",
    #     "pygraph",
    #     "subprocess",
    #     "dont_remove_include",
    # )

    dri = dont_remove_include.DontRemoveInclude()
    # sixth, iterate across all equivalence sets
    count_round = 0
    for es in equiv_sets:
        count_round += 1

        es_filtered = list(filter(dri.attempt_include_removal_for_file, es))
        random.shuffle(es_filtered)
        nfiles_to_process = len(es_filtered)
        nfiles_per_cpu = int(math.ceil(nfiles_to_process / fork_manager.max_n_jobs))
        print("Starting round", count_round, "with", nfiles_per_cpu, "jobs per cpu")
        sys.stdout.flush()
        es_subsets = []
        start = 0
        for i in range(fork_manager.max_n_jobs - 1):
            es_subsets.append(es_filtered[start : start + nfiles_per_cpu])
            start += nfiles_per_cpu
        es_subsets.append(es_filtered[start:])
        for count_jobid, es_subset in enumerate(es_subsets):
            pid = fork_manager.mfork()
            if pid == 0:
                trim_inclusions_from_files_extreme(es_subset+1, count_jobid)
                sys.exit(0)

        fork_manager.wait_for_remaining_jobs()
        tar_everything("bu_wholeshebang_round_" + str(count_round))
