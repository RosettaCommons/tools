#!/usr/bin/env python

######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################

from FARFAR_dag_util import *

######################################################################

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
import os
from os.path import basename, dirname, exists, expanduser
from time import sleep
from sets import Set
import copy

######################################################################

python_command=list_to_string(argv)
print_title_text("Enter: " + python_command)


DAG_ID=parse_options( argv, "DAG_ID", -1 )
num_struct_kept=parse_options(argv, "num_struct_kept", 0 )
filter_low_RMSD=parse_options(argv, "filter_low_RMSD", "True")

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )

delete_files=True

if(DAG_ID<0): error_exit_with_message("DAG_ID (%s) < 0! " %(DAG_ID))

if(num_struct_kept<=0): error_exit_with_message("num_struct_kept<=0! " %(num_struct_kept))

DAG_FOLDER=get_DAG_ID_folder(DAG_ID)

if(exists(DAG_FOLDER)==False): submit_subprocess("mkdir %s" %(DAG_FOLDER))

cat_outfile="%s/catfile.out" %(DAG_FOLDER)

######################################################################################
deleting_files_FARFAR_filtering_signal_file="%s/deleting_files_DAG_%d_FARFAR_filtering_signal.txt" %(DAG_FOLDER, DAG_ID)

if(exists(deleting_files_FARFAR_filtering_signal_file)):	error_exit_with_message("%s exist!" %(deleting_files_FARFAR_filtering_signal_file))

###########################################
#concatenate all the silent_files..

easy_cat_command='FARFAR_easy_cat.py -glob_folders %s -cat_file %s' %(DAG_FOLDER, cat_outfile ) 

print easy_cat_command
submit_subprocess(easy_cat_command)

############################################
filtered_RMSD_file=get_filtered_RMSD_file(DAG_ID)
filtered_energy_file=get_filtered_energy_file(DAG_ID)

if(exists(filtered_RMSD_file)): error_exit_with_message("filtered_RMSD_file (%s) already exist!" %(filtered_RMSD_file) )  
if(exists(filtered_energy_file)): error_exit_with_message("filtered_energy_file (%s) already exist!" %(filtered_energy_file) ) 

filter_energy_command="SWA_filter_outfile_wrapper.py -infile %s -scorecol_name score -max_n_struct %s -filter_outfile %s -remove_SCORE_file True " %(cat_outfile, num_struct_kept, filtered_energy_file)

print filter_energy_command
submit_subprocess(filter_energy_command)

if(filter_low_RMSD):
	filter_RMSD_command="SWA_filter_outfile_wrapper.py -infile %s -scorecol_name Full_L_rmsd -max_n_struct %s -filter_outfile %s -remove_SCORE_file True " %(cat_outfile, num_struct_kept, filtered_RMSD_file)

	print filter_RMSD_command
	submit_subprocess(filter_RMSD_command)

SCORE_silentfile='%s/%s' %(DAG_FOLDER, get_SCORE_file(DAG_ID))
create_SCORE_file_command = 'grep "SCORE: " %s > %s' %(cat_outfile, SCORE_silentfile)
print "create_SCORE_file_command= %s " %(create_SCORE_file_command) 
submit_subprocess(create_SCORE_file_command)
############################################

######################################################################################

if(delete_files): 

	create_generic_done_signal_file(deleting_files_FARFAR_filtering_signal_file)	#The mark the point where the scipt is now not rerunnable!

	print "rm -r %s" %(cat_outfile)
	submit_subprocess("rm -r %s" %(cat_outfile))

	###DON'T DELETE THE DAG_ID/ FOLDER, but instead delete the subfolders inside it!
	globstring =  '%s/*/silent_file.out' %(get_DAG_ID_folder(DAG_ID))

	globfiles = glob( globstring )
	glob_folders = map( lambda x : dirname(x) , globfiles )

	glob_folders.sort()
	for folder in glob_folders:
		command = 'rm -rf '+ folder
		print( command )
		sys.stdout.flush()
		sys.stderr.flush()
		submit_subprocess( command ) 
		
print_title_text("Exit: " + python_command)
