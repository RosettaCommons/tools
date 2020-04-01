import os
import sys

sys.path.insert( 0, os.path.realpath(__file__).rpartition("/")[0]+"/../python_cc_reader" )

from python_cc_reader.code_utilities import scan_compilable_files,compiled_cc_files


os.chdir("/Users/sheffler/hg/rosetta/rosetta_source/src")


files = set(scan_compilable_files())
files |= set(compiled_cc_files())
files = list(files)
files.sort()
for f in files:
	print f
