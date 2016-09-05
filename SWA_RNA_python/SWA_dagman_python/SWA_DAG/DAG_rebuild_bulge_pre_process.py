#!/usr/bin/env python


from os import system,popen
from os.path import exists,dirname,basename,expanduser,abspath
from sys import exit, argv
import string
from time import sleep
import glob
import time
import fileinput
import os
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
from SWA_dagman_python.utility.DAGMAN_util import update_CONDOR_file_with_actual_job_queue_num
######################################################################


###Example: /home/sripakpa/SWA_dagman_python/DAG_rebuild_bulge_pre_process.py -silent_file region_FINAL.out -condir_submit_file CONDOR/REBUILD_BULGE_FINAL.condor 


#run by master node, but through a submit_process()

python_command = get_python_command(argv)
print_title_text("Enter: " + python_command)

input_silent_file = parse_options( argv, "silent_file", "" ) 

if(input_silent_file==""): error_exit_with_message("User need to pass in input_silent_file ")

if(exists( input_silent_file )==False): error_exit_with_message("input_silent_file (%s) doesn't exist!" %(abspath(input_silent_file) ))

condor_submit_file = parse_options( argv, "condor_submit_file", "" )  

if(condor_submit_file==""): error_exit_with_message("User need to pass in condor_submit_file")

if(exists( condor_submit_file )==False): error_exit_with_message("condor_submit_file (%s) doesn't exist!" %(abspath(condor_submit_file) ))

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )


######TO DO:CHECK if the input silent_file is not empty###############


#################################################################################################################
data=safe_open(input_silent_file, mode='r', Is_master=False)

SEQUENCE_LINE=data.readline()
COLUMN_NAME_LINE=data.readline()

COL_NAME_LIST=COLUMN_NAME_LINE.split()

try:
	tag_col_index=COL_NAME_LIST.index('description')
except:
	print "COL_NAME_LIST=", COL_NAME_LIST
	error_exit_with_message("Cannot find description column index!")

if(tag_col_index!=(len(COL_NAME_LIST)-1)): error_exit_with_message("tag_col_index!=(len(COL_NAME_LIST)-1)")

print "tag is located at column_num: %d" %(tag_col_index+1)

pose_num=0

while(True):

	line=data.readline()

	if(line==''): break #End of file!

	if(len(line) <= 1): error_exit_with_message("len(line) <= 1") #check for line with only '\n'

	if( (line.find("SCORE:") != -1) != (line[0:6]=="SCORE:") ): error_exit_with_message("(line.find(\"SCORE:\") != -1) != (line[0:6]==\"SCORE:\") for line=%s" %(line))

	if(line[0:6]=="SCORE:"): #SCORE line! 
		
		if(line.find('description') != -1): error_exit_with_message("extra column_name line (%s)" %(line) )
	
		pose_num+=1
	
		score_line=line

		cols=score_line.split()

		actual_tag=cols[tag_col_index]

		assumed_tag = 'S_%d' %(pose_num-1) #Oct 22, 2011: Assume this form for now

		if(assumed_tag!=actual_tag): error_exit_with_message("assumed_tag (%s)!=actual_tag (%s)" %(assumed_tag, actual_tag) )	

data.close()

N_JOBS = pose_num

update_CONDOR_file_with_actual_job_queue_num(condor_submit_file, N_JOBS)

print_title_text("Exit: " + python_command)



