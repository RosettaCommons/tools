#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
import os
import copy

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options



input_sequence = parse_options( argv, "input_sequence", "" )

lower_sequence = input_sequence.lower()

print "lower_sequence=%s " %(lower_sequence)
