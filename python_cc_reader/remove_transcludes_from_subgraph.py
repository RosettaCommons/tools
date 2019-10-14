from python_cc_reader.inclusion_removal.remove_headers_from_subgraph import remove_transcludes_from_subgraph
from python_cc_reader.utility.fork_manager import ForkManager
from python_cc_reader.external.blargs import blargs

if __name__ == "__main__":
    with blargs.Parser(locals()) as p:
        p.int("ncpu")
        p.require_one(
            p.multiword("starting_files").cast(lambda x: x.split()),
            p.str("starting_file_list")
        )

    # 1. Load the source tree
    # 2. Construct the transitive closure graph
    # 3. Add transitive includes to file X that are
    #    included or transcluded by file Y
    #    listed in the starting_files set that X
    #    includes or transcludes.
    # 4. Construct the equivalence sets
    # 5. In batches, process the file X that transclude one
    #    of the starting files
    #    -- try to remove the headers that were transcluded
    #       by the header Y in the original inclusion graph
    #       that X transcludes


    fork_manager = ForkManager(ncpu)

    if starting_file_list is not None:
        starting_files = [line.strip() for line in open(starting_file_list).readlines()]

    remove_transcludes_from_subgraph(starting_files, fork_manager)

