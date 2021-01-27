from python_cc_reader.inclusion_removal.remove_headers_from_all_files import create_inclusion_graph_from_source
from python_cc_reader.inclusion_removal.inclusion_graph import (
    transitive_closure,
    add_indirect_headers,
    trim_unnecessary_headers_from_file,
    wrap_cxxtest_compile,
    canonical_prohibit_removal,
    total_order_from_graph,
)
from python_cc_reader.cpp_parser.code_utilities import compiled_cxxtest_hh_files
from python_cc_reader.inclusion_removal.add_headers import write_file
from python_cc_reader.inclusion_removal import remove_headers_from_all_files
from python_cc_reader.inclusion_removal.dont_remove_include import DontRemoveInclude
from python_cc_reader.external.blargs import blargs
from python_cc_reader.utility.fork_manager import ForkManager
import sys

if __name__ == "__main__":
    with blargs.Parser(locals()) as p:
        p.int("ncpu").shorthand("j").default(1)
    
    fork_manager = ForkManager(ncpu)
    dri = DontRemoveInclude()
    fnames = compiled_cxxtest_hh_files()

    compilable_files, all_includes, file_contents, g = (
        create_inclusion_graph_from_source("pickled_include_graph.bin")
    )
    for cxxfile in fnames:
        file_contents[cxxfile] = open(cxxfile).readlines()

    tg = transitive_closure(g)
    total_order = total_order_from_graph(g)

    print("beginning examination of", len(fnames), "cxxtest files")

    for i, fname in enumerate(fnames):
       pid = fork_manager.mfork()
       if pid == 0:
          id = fork_manager.myjob_index+1
          trim_unnecessary_headers_from_file(
             fname,
             tg,
             total_order,
             file_contents,
             wrap_cxxtest_compile(id),
             canonical_prohibit_removal(dri),
          )
          write_file(fname, file_contents[fname])
          sys.exit(0)
    fork_manager.wait_for_remaining_jobs()

