import sys

from python_cc_reader.inclusion_removal import remove_headers_from_all_files
from python_cc_reader.external.blargs import blargs
from python_cc_reader.utility.fork_manager import ForkManager

if __name__ == "__main__":
    with blargs.Parser(locals()) as p:
        p.int("ncpu").shorthand("j").default(1)
        p.flag("skip_1st_compile")
        p.int("starting_round").default(1)
    
    fork_manager = ForkManager(ncpu)
    remove_headers_from_all_files.whole_shebang(
        fork_manager, skip_1st_compile, starting_round)
