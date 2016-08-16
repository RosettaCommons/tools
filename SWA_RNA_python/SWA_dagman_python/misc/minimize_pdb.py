#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options


#minimize_pdb.py -pdb_list XX -minimize_res XX

pdb_list =	  parse_options( argv, "pdb_list", [""] )
minimize_res =	  parse_options( argv, "minimize_res", [ 0 ] )
force_field_file = parse_options( argv, "force_field_file", "")
fold_tree_strings = parse_options( argv, "fold_tree_strings", [""])

if(force_field_file==""): error_exit_with_message("User need to pass in force_field_file!")
if(len(minimize_res)==0): error_exit_with_message("User need to pass in minimize_res!")
if(pdb_list==[""]): error_exit_with_message("User need to pass in moving_pdb_list!")

if(len(argv)!=1):
	print argv, " leftover len(argv)=", len(argv)
	assert(False)

if(fold_tree_strings!=[""]):
	for n in range(len(fold_tree_strings)):
		fold_tree_strings[n]=fold_tree_strings[n].replace("_","-")

##############################################################
HOMEDIR = expanduser('~') 

EXE = HOMEDIR+'/src/mini/bin/parin_test.linuxgccrelease'

if not( exists( EXE )):
	EXE = HOMEDIR+'/src/mini/bin/parin_test.macosgccrelease'

print "HOMEDIR= %s" %(HOMEDIR)
assert( exists( EXE ) )

DB = HOMEDIR+'/minirosetta_database'
assert( exists( DB ) )
##############################################################


command=EXE
command += ' -algorithm minimize_pdb'
command += ' -database ~/minirosetta_database'
command += ' -s ' + list_to_string(pdb_list)
command += ' -minimize_res %s' %(list_to_string(minimize_res))
command += ' -force_field_file %s' %(force_field_file)

if(fold_tree_strings!=[""]): command += ' -fold_tree_strings %s ' %(list_to_string(fold_tree_strings))

command += ' > minimize_pdb_output.txt '

print command 
submit_subprocess( command )

