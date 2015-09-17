# This is a script meant to be run from the command line
# that allows a user to list of .cc files that the user
# wants to trim unnecessary #includes from.
# The script, just like all the python code utilities
# is meant to be run from within the mini/src/ directory
# usage:
#
# python trim_unnecessary_headers_from_cc.py <cc_filename> [<cc_filename>]*



from .inclusion_graph import *
import sys

if len( sys.argv ) < 2:
    print("Error: expected a .cc file list to remove headers from")
    sys.exit(1)

flist = sys.argv[ 1: ]
for fname in flist :
    print("Adding", fname, "to the list of files to trim")
    if len( fname ) < 3 :
        print("Error: the filename", fname,"is too short!")
        sys.exit(1)
    if fname[:-3] == ".hh" or (len(fname)>3 and fname[:-4] == ".hpp" ) :
        print("ERROR: this script is only useful for trimming .cc files.")
        print("File", fname, "cannot be processed")
        sys.exit(1)

trim_inclusions_from_files_extreme(flist, 1)
