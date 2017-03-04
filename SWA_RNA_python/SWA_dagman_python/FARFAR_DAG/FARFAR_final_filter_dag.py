#!/usr/bin/env python

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

from FARFAR_dag_util import *
######################################################################

from SWA_dagman_python.misc.SWA_cat_outfiles import *
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################




python_command=list_to_string(argv)
print_title_text("Enter: " + python_command)

num_DAG= parse_options( argv, "num_DAG", 0 )
filter_low_RMSD=parse_options(argv, "filter_low_RMSD", "True")

delete_files=False

if(num_DAG <= 0 ): error_exit_with_message('user need to pass in a positive integer for num_DAG')

final_SCORE_file="FINAL_SCORE.out"

final_filtered_RMSD_file=parse_options(argv, "final_filtered_RMSD_file", "")
final_filtered_energy_file=parse_options(argv, "final_filtered_energy_file", "")

if(final_filtered_RMSD_file==""): error_exit_with_message("final_filtered_RMSD_file==\"\"")
if(final_filtered_energy_file==""): error_exit_with_message("final_filtered_energy_file==\"\"")

if(exists(final_SCORE_file)): error_exit_with_message('final_SCORE_file (%s) already exist!' %(final_SCORE_file))
if(exists(final_filtered_RMSD_file)): error_exit_with_message('final_filtered_RMSD_file (%s) already exist!' %(final_filtered_RMSD_file))
if(exists(final_filtered_energy_file)): error_exit_with_message('final_filtered_energy_file (%s) already exist!' %(final_filtered_energy_file))


filtered_RMSD_file_list=[]
filtered_energy_file_list=[]
SCORE_file_list=[]

for DAG_ID in range(num_DAG):
	
	DAG_FOLDER=get_DAG_ID_folder(DAG_ID)

	if(exists(DAG_FOLDER)==False): error_exit_with_message("DAG_FOLDER (%s) doesn't exist! " %(DAG_FOLDER))

	#combine the score file
	SCORE_file='%s/%s' %(DAG_FOLDER, get_SCORE_file(DAG_ID) )

	if(exists(SCORE_file)==False): error_exit_with_message("SCORE_file (%s) doesn't exist! " %(SCORE_file))

	SCORE_file_list.append(SCORE_file)

	submit_subprocess("cat %s >> %s " %(SCORE_file, final_SCORE_file)) #does this have a line limit?

	filtered_energy_file=get_filtered_energy_file(DAG_ID)

	if(exists(filtered_energy_file)==False): error_exit_with_message("filtered_energy_file (%s) doesn't exist! " %(filtered_energy_file))

	filtered_energy_file_list.append(filtered_energy_file)

	if(filter_low_RMSD):

		filtered_RMSD_file=get_filtered_RMSD_file(DAG_ID)

		if(exists(filtered_RMSD_file)==False): error_exit_with_message("filtered_RMSD_file (%s) doesn't exist! " %(filtered_RMSD_file))

		filtered_RMSD_file_list.append(filtered_RMSD_file)


print "filtered_energy_file_list: ", filtered_energy_file_list
concatenate_outfiles(infile_list=filtered_energy_file_list, outfile=final_filtered_energy_file)

if(filter_low_RMSD):
	print "filtered_RMSD_file_list: ", filtered_RMSD_file_list
	concatenate_outfiles(infile_list=filtered_RMSD_file_list, outfile=final_filtered_RMSD_file)


if(delete_files):
	for filtered_energy_file in filtered_energy_file_list:
		submit_subprocess("rm %s " %(filtered_energy_file)) #does this have a line limit?

	for filtered_RMSD_file in filtered_RMSD_file_list:
		submit_subprocess("rm %s " %(filtered_RMSD_file)) #does this have a line limit?


print_title_text("Exit: " + python_command)


#sripakpa@sripakpas-MacBook-Pro:~/minirosetta/11_2010/Nov_4_FARFAR_simple_motif_allow_bulge_mode/main_folder/1S72_GAAA$ du -sh SCORE_cat_file.out 
# 35M    SCORE_cat_file.out

#sripakpa@sripakpas-MacBook-Pro:~/minirosetta/11_2010/Nov_4_FARFAR_simple_motif_allow_bulge_mode/main_folder/1S72_GAAA$ du -sh cat_file.out 
#439M    cat_file.out

