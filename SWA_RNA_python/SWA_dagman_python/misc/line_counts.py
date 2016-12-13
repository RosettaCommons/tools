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


filename=parse_options( argv, "file", "")

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_args=%s " %(list_to_string(argv) ) )

if(filename==""): error_exit_with_message('User need to specify -file option!')

if(exists(filename)==False): error_exit_with_message("filename (%s) doesn't exist!" %(filename) )
	

FILE = open( filename, 'r')

line_count=0
word_count=0
char_count=0

for line in FILE:

	line_count+=1
	word_count+=len(line.split())
	char_count+=len(line)


FILE.close()

print "line_count=%s | word_count=%s | char_count=%s" %(line_count, word_count, char_count)
