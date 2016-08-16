#!/usr/bin/env python

######################################################################

from SWA_dagman_python.utility.DAGMAN_util import *
from SWA_dagman_python.utility.SWA_util import *

from SWA_dagman_python.parser.SWA_parse_options import replace_arg_value

######################################################################

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
import os
import copy

from sets import Set


######################################################################

DAG_FILE_FOLDER="CONDOR/" 

MINIMIZER_FOLDER="FULL_LENGTH_MINIMIZER"

NSTRUCT_PER_MAPPER=10

GENERIC_SILENT_FILE_REDUCER_SCRIPT= get_PYEXE("dagman/DAG_generic_silent_file_reducer.py")

SWA_RNA_MAIN_EXE=""
if(use_new_src_code()):
	SWA_RNA_MAIN_EXE= get_rosetta_EXE("swa_rna_main") 
	error_exit_with_message("rna_minimize algorithm not yet implemented in swa_rna_main!");
else:
	SWA_RNA_MAIN_EXE= get_rosetta_EXE("rna_swa_test") 


######################################################################
def get_rna_minimize_process_output_filename():
	return "%s/pre_process_tag_num_range.txt" %(MINIMIZER_FOLDER)



######################################################################
def swa_minimize_preprocess(in_silent_file, in_tag_prefix, tag_range, minimizer_dag_job_file):

	###################################preprocess DAG_FILE################################
	
	data=safe_open(in_silent_file, mode='r', Is_master=False)

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



	tag_info_list = []

	while(True):

		line=data.readline()

		if(line==''): break #End of file!

		if(len(line) <= 1): error_exit_with_message("len(line) <= 1") #check for line with only '\n'

		if( (line.find("SCORE:") != -1) != (line[0:6]=="SCORE:") ): error_exit_with_message("(line.find(\"SCORE:\") != -1) != (line[0:6]==\"SCORE:\") for line=%s" %(line))

		if(line[0:6]=="SCORE:"): #SCORE line! 
		
			if(line.find('description') != -1): error_exit_with_message("extra column_name line (%s)" %(line) )
	
	
			score_line=line

			cols=score_line.split()

			actual_tag=cols[tag_col_index]

			tag_info={}

			tag_info["ID"]=int(actual_tag.split("_")[-1])

			assumed_tag="%s%d" %(in_tag_prefix, tag_info["ID"])

			if(assumed_tag!=actual_tag): error_exit_with_message("assumed_tag=(%s)!=(%s)=actual_tag" %(assumed_tag, actual_tag) )

			tag_info["tag"]=actual_tag

			tag_info_list.append(tag_info)


	data.close()


	########################################################################################
	min_tag_num=0
	max_tag_num=0

	if(len(tag_range)!=0):	

		if(len(tag_range)!=2): error_exit_with_message("tag_range (%s) is not empty but len(tag_range)!=2" %( list_to_string(tag_range) ) )
		min_tag_num=tag_range[0]
		max_tag_num=tag_range[1]
		if(min_tag_num>max_tag_num): error_exit_with_message("min_tag_num>max_tag_num")

		print "min_tag_num=%s | max_tag_num=%d " %(min_tag_num, max_tag_num)

	########################################################################################

	tag_info_list=sort_dictionary_list(tag_info_list, key="ID") 
	filtered_tag_info_list=[]
	
	for silent_struct_ID in range(len(tag_info_list)):

		tag_info=tag_info_list[silent_struct_ID]

		if(len(tag_range)!=0):	
			if(tag_info["ID"]<min_tag_num): continue
			if(tag_info["ID"]>max_tag_num): continue

		filtered_tag_info_list.append(copy.deepcopy(tag_info))

	tag_info_list=copy.deepcopy(filtered_tag_info_list)

	########################################################################################
	pre_process_output_filename=get_rna_minimize_process_output_filename() 

	if(exists(pre_process_output_filename)): 
		print "pre_process_output_filename (%s) already exist! ... removing..." %(pre_process_output_filename)

	outfile = safe_open(pre_process_output_filename, mode='w' ,Is_master=False)

	########################################################################################
	curr_struct_count=0
	queue_ID=-1

	for silent_struct_ID in range(len(tag_info_list)):

		if(curr_struct_count==0): 
			queue_ID+=1
			outfile.write("%10s " %(queue_ID))

		tag_info=tag_info_list[silent_struct_ID]

		curr_struct_count+=1
	
		outfile.write("%10s " %(tag_info["tag"]))

		if( (curr_struct_count==NSTRUCT_PER_MAPPER) or (silent_struct_ID==(len(tag_info_list)-1)) ):
			outfile.write(" \n")
			curr_struct_count=0

			
	N_JOBS=queue_ID

	print "N_JOBS=%s | NSTRUCT_PER_MAPPER=%d | len(tag_info_list)=%d" %(N_JOBS, NSTRUCT_PER_MAPPER, len(tag_info_list))

	outfile.close()
	########################################################################################


	update_CONDOR_file_with_actual_job_queue_num(minimizer_dag_job_file, N_JOBS)

	########################################################################################
