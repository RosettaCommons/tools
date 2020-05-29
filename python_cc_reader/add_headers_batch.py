#!/usr/bin/env python

import sys

# sys.path.append("/Users/brian/SVN/mini_tools")
from python_cc_reader.inclusion_removal import add_headers

from optparse import OptionParser

usage = "%prog file_list header_name"
parser = OptionParser()
(options, args) = parser.parse_args()

file = open(args[0], "r")

file_list = []
for line in file:
    file_list.append(line.rstrip())

for path in file_list:
    print(path)
    add_headers.add_autoheaders_to_file(path, [args[1]])
