#!/usr/bin/env python


import sys
import string
from os import system,popen
from os.path import basename, dirname, exists, expanduser
from glob import glob
import os
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
from SWA_dagman_python.misc.SWA_cat_outfiles import *
######################################################################

cat_file=parse_options( argv, "cat_file", "cat_file.out" )
#glob_pattern=parse_options( argv, "glob_pattern", ".out" ) #Commented this on Jan 14, 2011

glob_folders=parse_options( argv, "glob_folders", [""] )

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )

if(exists(cat_file)): 
	print "cat_file (%s) already exist... removing..." %(cat_file)
	print "cat_file (%s) already exist... removing..." %(os.path.abspath(cat_file))
	submit_subprocess('rm -r %s ' %(cat_file))

if(glob_folders==[""]): error_exit_with_message('glob_folders==[""]')

print "glob_folders= ", glob_folders

#outfiles = sys.argv[1:]

silent_file_list=[]

for glob_folder in glob_folders:

	globfiles = glob( glob_folder + '/*/silent_file.out' )

	for file_name in globfiles:
		silent_file_list.append( file_name )             						


print "silent_file_list= " , silent_file_list


print "------------------------------------------------------------------------------------------"
print "------------------------------------------------------------------------------------------"

print "Catting into: ", os.path.abspath(cat_file) 

concatenate_outfiles(infile_list=silent_file_list, outfile=cat_file) #change from FARFAR_cat_outfiles.py to SWA_cat_outfiles on Nov 4, 2010

if(len(silent_file_list)!=0):
	lines = popen( 'grep "SCORE: " '+cat_file).readlines()
	print '... from %d primary files. Found %d  decoys.' % (len( silent_file_list ),len(lines)-1)
else:
	print " silent_file_list is still empty! "

print "------------------------------------------------------------------------------------------"
print "------------------------------------------------------------------------------------------"

