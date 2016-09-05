#!/usr/bin/env python

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
import os
from os.path import basename, dirname, exists, expanduser
from time import sleep
import copy
import string
from os import popen 
import math

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
from SWA_dagman_python.utility.RNA_sequence import get_one_letter_code


#reweight_chem_shift_score.py -infile minimize_WITH_SHIFT_STATS_Feb_06_SWA.out -new_chem_shift_weight 2.0 -use_rosetta_columns true



input_silent_file=parse_options( argv, "infile", "" )

if(input_silent_file==""): error_exit_with_message('input_silent_file==""')

if(exists(input_silent_file)==False): error_exit_with_message("input_silent_file (%s) doesn't exist!" %(input_silent_file))

new_chem_shift_weight=parse_options( argv, "new_chem_shift_weight", -999999.99)

if(new_chem_shift_weight<-999998.88): error_exit_with_message("User need to pass in new_chem_shift_weight option!")

output_silent_file=parse_options( argv, "outfile", "" )



use_rosetta_columns=parse_options( argv, "use_rosetta_columns", "")

num_shift_data_col_name=""
shift_RMSD_col_name=""
shift_score_col_name=""

new_outfile_chem_shift_str=""

if(use_rosetta_columns=="true"):
	num_shift_data_col_name="R_num_shift_data"
	shift_RMSD_col_name="R_shift_RMSD"
	shift_score_col_name="rna_chem_shift"
	new_outfile_chem_shift_str="NEW_RCHEM_SHIFT"
elif(use_rosetta_columns=="false"):
	num_shift_data_col_name="num_shift_data"
	shift_RMSD_col_name="shift_RMSD"
	shift_score_col_name="shift_score"
	new_outfile_chem_shift_str="NEW_CHEM_SHIFT"
else:
	error_exit_with_message("Invalid use_rosetta_columns (%s) option!" %(use_rosetta_columns))

#####################################################################

if(output_silent_file==""):
	output_silent_file="%.1fX_%s_%s" %(new_chem_shift_weight, new_outfile_chem_shift_str, basename(input_silent_file))

else:
	if(exists(output_silent_file)):
		error_exit_with_message("User-specified outfile (%s) already exist!" %(output_silent_file))

print "output_silent_file=%s" %(output_silent_file)

#####################################################################

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )
#####################################################################
data=safe_open(input_silent_file, mode='r', Is_master=False)

SEQUENCE_LINE=data.readline()
COLUMN_NAME_LINE=data.readline()
COL_NAME_LIST=COLUMN_NAME_LINE.split()

#####################################################################

try:
	num_shift_data_index=COL_NAME_LIST.index(num_shift_data_col_name)
except:
	print "COL_NAME_LIST=", COL_NAME_LIST
	error_exit_with_message("Cannot find column index of num_shift_data_col_name (%s)!" %(num_shift_data_col_name) )


try:
	shift_RMSD_index=COL_NAME_LIST.index(shift_RMSD_col_name)
except:
	print "COL_NAME_LIST=", COL_NAME_LIST
	error_exit_with_message("Cannot find column index of shift_RMSD_col_name (%s)!" %(shift_RMSD_col_name) )

try:
	shift_score_index=COL_NAME_LIST.index(shift_score_col_name)
except:
	print "COL_NAME_LIST=", COL_NAME_LIST
	error_exit_with_message("Cannot find column index of shift_score_col_name (%s)!" %(shift_score_col_name) )


try:
	score_col_index=COL_NAME_LIST.index('score')
except:
	print "COL_NAME_LIST=", COL_NAME_LIST
	error_exit_with_message("Cannot find score column index!" )

if(score_col_index!=1): error_exit_with_message("score_col_index!=1")

try:
	tag_col_index=COL_NAME_LIST.index('description')
except:
	print "COL_NAME_LIST=", COL_NAME_LIST
	error_exit_with_message("Cannot find description column index!")

if(tag_col_index!=(len(COL_NAME_LIST)-1)): error_exit_with_message("tag_col_index!=(len(COL_NAME_LIST)-1)")


assert_no_duplicate_in_string_list(COL_NAME_LIST)

#####################################Consistency check################################################
if(COL_NAME_LIST.count("SCORE:")!=1): error_exit_with_message("COL_NAME_LIST.count(\"SCORE:\")!=1")
if(COL_NAME_LIST.index("SCORE:")!=0): error_exit_with_message("COL_NAME_LIST.index(\"SCORE:\")!=0")

if(COL_NAME_LIST.count("description")!=1): error_exit_with_message("COL_NAME_LIST.count(\"description\")!=1")
if(COL_NAME_LIST.index("description")!=(len(COL_NAME_LIST)-1)): 
	error_exit_with_message("COL_NAME_LIST.index(\"description\")!=(len(COL_NAME_LIST)-1)")

############4. Create the new silent_file with new chem_shift_score#######



OUTPUT_SILENT_FILE=open(output_silent_file, 'w')
 
OUTPUT_SILENT_FILE.write(SEQUENCE_LINE) #Sequence line

OUTPUT_SILENT_FILE.write(COLUMN_NAME_LINE)

######################################################################################################

######################################################################################################
first_update=True

while(True):

	line=data.readline()

	if(line==''): break #End of file!

	if(len(line) <= 1): error_exit_with_message("len(line) <= 1") #check for line with only '\n'

	if( (line.find("SCORE:") != -1) != (line[0:6]=="SCORE:") ): error_exit_with_message("(line.find(\"SCORE\") != -1) != (line[0:6]==\"SCORE:\")")

	if(line[0:6]=="SCORE:"): #SCORE line! 
		
		if(line.find('description') != -1): error_exit_with_message("extra column_name line (%s)" %(line) )
		
		score_line=line

		cols=score_line.split()

		total_score=float(cols[score_col_index])

		shift_RMSD=float(cols[shift_RMSD_index])
		num_shift_data=float(cols[num_shift_data_index])

		old_shift_score=float(cols[shift_score_index])
		new_shift_score=new_chem_shift_weight*(shift_RMSD*shift_RMSD)*num_shift_data

		if(first_update):
			first_update=False
			print "OLD_shift_weight is %.1f" %(old_shift_score/((shift_RMSD*shift_RMSD)*num_shift_data))

		cols[shift_score_index]="%8.3f" %( new_shift_score ) #Dec 10, 2011: USED TO BE "%10.3f" %(total_score) #OK since include extra white_space character below

		total_score=total_score-old_shift_score+new_shift_score

		cols[score_col_index]="%8.3f" %(total_score) #Dec 10, 2011:"%10.3f" %(total_score) #OK since include extra white_space character below

		if(len(cols)!=len(COL_NAME_LIST)): error_exit_with_message("len(cols)!=len(COL_NAME_LIST)")

		new_score_line=""		

		for n in range(len(COL_NAME_LIST)):
			if(n==0):
				new_score_line+="%-*s " %(len(COL_NAME_LIST[n])+2, cols[n])
			else:
				new_score_line+=" %*s"  %(max(len(COL_NAME_LIST[n])+2,10), cols[n])

		if(len(new_score_line.split())!=len(COL_NAME_LIST)): 
			error_exit_with_message("len(new_score_line.split())=(%s)!=(%s)=len(COL_NAME_LIST) for new_score_line=%s" %(len(new_score_line.split()),len(COL_NAME_LIST), new_score_line))

		OUTPUT_SILENT_FILE.write(new_score_line +'\n')

	else:
		OUTPUT_SILENT_FILE.write(line)

	offset=data.tell()

data.close()

OUTPUT_SILENT_FILE.close()


