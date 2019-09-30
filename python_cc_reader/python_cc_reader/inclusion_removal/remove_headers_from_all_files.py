import subprocess
import math
import random
import re
import sys
import pickle
import os

from ..cpp_parser.code_utilities import load_source_tree, expand_includes_for_file

from .inclusion_graph import create_graph_from_includes, remove_known_circular_dependencies_from_graph, transitive_closure, add_indirect_headers, inclusion_equivalence_sets_from_desired_subgraph, write_equiv_sets_file, trim_inclusions_from_file_extreme, total_order_from_graph
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


def verify_all_files_compile(fork_manager, compilable_files):
    d = {}
    d["any_fail_to_compile"] = False
    def set_any_failed(_1, _2):
        print("FAILURE CALLBACK");
        d["any_fail_to_compile"] = True
    fork_manager.error_callback = set_any_failed
    
    for fname in compilable_files :
        if d["any_fail_to_compile"]:
            break
        pid = fork_manager.mfork()
        if pid == 0:
            if d["any_fail_to_compile"]:
                sys.exit(0)
            print("testing compilation of", fname)
            if not test_compile(fname) :
                print( "Error: ", fname, "does not compile on its own")
                sys.exit(1)
            else:
                sys.exit(0)
    fork_manager.wait_for_remaining_jobs()
    return d["any_fail_to_compile"]

def whole_shebang(fork_manager: ForkManager, skip_test_compile=False, starting_round=1):

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

    if os.path.isfile("pickled_include_graph.bin"):
        with open("pickled_include_graph.bin","rb") as fid:
            graphs = pickle.load(fid)
        compilable_files = graphs["compilable_files"]
        all_includes = graphs["all_includes"]
        file_contents = graphs["file_contents"]
        g = graphs["g"]
    else :
        # first create the include graph
        # trim this graph to remove known cycles (vector1_bool, vector0_bool)

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
        with open("pickled_include_graph.bin","wb") as fid:
            pickle.dump(graphs, fid)

    
    # second, verify that all files compile on their own

    if not skip_test_compile:
        any_fail_to_compile = verify_all_files_compile(fork_manager, compilable_files)
        if any_fail_to_compile:
            print("Error: coud not compile all files on their own")
            sys.exit(0)

    if os.path.isfile("pickled_equiv_sets.bin"):
        with open("pickled_equiv_sets.bin","rb") as fid:
            eqsets = pickle.load(fid)
        tg = eqsets["tg"]
        compilable_files = eqsets["compilable_files"]
        equiv_sets = eqsets["equiv_sets"]
    else :
        # third, compute the transitive closure graph
        print("computing transitive closure")
        tg = transitive_closure(g)
    
        # fourth, add all headers in transitive closure
        tar_everything("bu_starting_code")
        print("adding indirect headers")
        add_indirect_headers(tg, compilable_files)
        tar_everything("bu_transclose_headers_round_0")
    
        # fifth, compute the equivalence sets, focusing only on files in
        # core/ protocols/ devel/ and apps/
        print("computing equivalence sets from inclusion graph")
        equiv_sets = inclusion_equivalence_sets_from_desired_subgraph(g)
        write_equiv_sets_file("equivalence_sets_wholeshebang.txt", equiv_sets)

        eqsets = {}
        eqsets["tg"] = tg
        eqsets["compilable_files"] = compilable_files
        eqsets["equiv_sets"] = equiv_sets
        with open("pickled_equiv_sets.bin","wb") as fid:
            pickle.dump(eqsets, fid)


    dri = DontRemoveInclude()
    # sixth, iterate across all equivalence sets
    # for each equivalence set, fork ncpu processes
    # and each process will remove #includes from a
    # subset of the files in that equivalence set
    count_round = 0
    for es in equiv_sets:
        count_round += 1
        if count_round < starting_round :
            continue
        es_filtered = list(filter(dri.attempt_include_removal_for_file, es))
        nfiles_to_process = len(es_filtered)
        print("Starting round", count_round, "with", nfiles_to_process, "to process")
        sys.stdout.flush()

        if count_round > 1:
            print("loading source tree")
            compilable_files, all_includes, file_contents = load_source_tree()
            print("creating graph from includes")
            g = create_graph_from_includes(all_includes)
            remove_known_circular_dependencies_from_graph(g)
            print("creating transitive closure")
            tg = transitive_closure(g)
        
        print("computing total order")
        sys.stdout.flush()
        total_order = total_order_from_graph(g)

        # es_subsets = []
        # start = 0
        # for i in range(fork_manager.max_n_jobs - 1):
        #     es_subsets.append(es_filtered[start : start + nfiles_per_cpu])
        #     start += nfiles_per_cpu
        # es_subsets.append(es_filtered[start:])
        for i, fname in enumerate(es_filtered):
            pid = fork_manager.mfork()
            # pid = 0 # TEMP !
            if pid == 0:
                trim_inclusions_from_file_extreme(
                    fname,
                    (i+1) % fork_manager.max_n_jobs,
                    compilable_files,
                    all_includes,
                    file_contents,
                    g,
                    tg,
                    total_order
                )
                sys.exit(0)

        fork_manager.wait_for_remaining_jobs()
        tar_everything("bu_wholeshebang_round_" + str(count_round))
