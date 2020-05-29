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


def create_inclusion_graph_from_source(pickle_file_name=None):
    if pickle_file_name and os.path.isfile(pickle_file_name):
        with open(pickle_file_name,"rb") as fid:
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

        if pickle_file_name:
            graphs = {}
            graphs["compilable_files"] = compilable_files
            graphs["all_includes"] = all_includes
            graphs["file_contents"] = file_contents
            graphs["g"] = g
            with open(pickle_file_name,"wb") as fid:
                pickle.dump(graphs, fid)

    return compilable_files, all_includes, file_contents, g


def add_transitive_closure_headers(g, files_needing_transcludes, pickled_tg_and_eqset_fname=None):
    if pickled_tg_and_eqset_fname and os.path.isfile(pickled_tg_and_eqset_fname):
        with open(pickled_tg_and_eqset_fname,"rb") as fid:
            eqsets = pickle.load(fid)
        tg = eqsets["tg"]
        equiv_sets = eqsets["equiv_sets"]
    else :
        # third, compute the transitive closure graph
        print("computing transitive closure")
        tg = transitive_closure(g)

        # fourth, add all headers in transitive closure
        tar_everything("bu_starting_code")
        print("adding indirect headers")
        add_indirect_headers(tg, files_needing_transcludes)

        # fifth, compute the equivalence sets, focusing only on files in
        # core/ protocols/ devel/ and apps/
        print("computing equivalence sets from inclusion graph")
        equiv_sets = inclusion_equivalence_sets_from_desired_subgraph(g)
        write_equiv_sets_file("equivalence_sets_wholeshebang.txt", equiv_sets)

        if picked_tg_and_eqset_fname:
            eqsets = {}
            eqsets["tg"] = tg
            eqsets["equiv_sets"] = equiv_sets
            with open(pickled_tg_and_eqset_fname,"wb") as fid:
                pickle.dump(eqsets, fid)

    return tg, equiv_sets

def remove_includes_from_files_in_parallel(
        fnames,
        fork_manager,
        compilable_files,
        all_includes,
        file_contents,
        g,
        tg,
        total_order,
        dri
):
    for i, fname in enumerate(fnames):
        pid = fork_manager.mfork()
        if pid == 0:
            trim_inclusions_from_file_extreme(
                fname,
                fork_manager.myjob_index+1,
                file_contents,
                tg,
                total_order,
                dri
            )
            sys.exit(0)
    fork_manager.wait_for_remaining_jobs()



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

    # first, load the inclusion graph
    compilable_files, all_includes, file_contents, g = (
        create_inclusion_graph_from_source("pickled_include_graph.bin")
    )

    # second, verify that all files compile on their own
    if not skip_test_compile:
        any_fail_to_compile = verify_all_files_compile(fork_manager, compilable_files)
        if any_fail_to_compile:
            print("Error: coud not compile all files on their own")
            sys.exit(0)

    # third, create transitive closure graph and determine file equivalence
    # sets from the partial order
    tg, equiv_sets = add_transitive_closure_headers(
        g, compilable_files, "pickled_equiv_sets.bin")
    tar_everything("bu_transclose_headers_round_0")

    dri = DontRemoveInclude()
    # sixth, iterate across all equivalence sets
    # for each equivalence set, fork ncpu processes
    # and each process will remove #includes from a
    # subset of the files in that equivalence set

    for count_round, es in enumerate(equiv_sets):
        if count_round+1 < starting_round :
            continue
        equiv_set_filtered = list(filter(dri.attempt_include_removal_for_file, es))
        nfiles_to_process = len(equiv_set_filtered)
        print("Starting round", count_round+1, "with", nfiles_to_process, "to process")
        sys.stdout.flush()

        if count_round > 1:
            compilable_files, all_includes, file_contents, g = (
                create_inclusion_graph_from_source()
            )
            print("creating transitive closure")
            tg = transitive_closure(g)
            print("number of edges in transitive closure graph:", len(tg.edge_properties))
        else:
            # the file contents have changed after we added the transcludes
            print("reloding file contents")
            file_contents = {fname: open(fname).readlines() for fname in file_contents.keys()}

        print("computing total order")
        sys.stdout.flush()
        total_order = total_order_from_graph(g)

        remove_includes_from_files_in_parallel(
            equiv_set_filtered,
            fork_manager,
            compilable_files,
            all_includes,
            file_contents,
            g,
            tg,
            total_order,
            dri
        )

        tar_everything("bu_wholeshebang_round_" + str(count_round))
