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


#reweight_silent_file_SCORE.py -infile WITH_SHIFT_STATS_Oct_20_SWA_Jan_28_1JJ2.out -reweight_list shift_score-0.0

#reweight_silent_file_SCORE.py -infile WITH_SHIFT_STATS_combine.out -reweight_list shift_score-3.0 

#reweight_silent_file_SCORE.py -infile WITH_SHIFT_STATS_March_04_SWA.out -reweight_list shift_score-1.0/3.0

#reweight_silent_file_SCORE.py -infile region_FINAL.out -reweight_list fa_stack-0.25  rna_bulge-4.5/10.0


input_silent_file=parse_options( argv, "infile", "" )

if(input_silent_file==""): error_exit_with_message('input_silent_file==""')

if(exists(input_silent_file)==False): error_exit_with_message("input_silent_file (%s) doesn't exist!" %(input_silent_file))

reweight_string_list=parse_options( argv, "reweight_list", [""]) #[column_name-weight_factor]..example:  fa_stack-4.0  rna_bulge-10.0/4.5

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )
#####################################################################
data=safe_open(input_silent_file, mode='r', Is_master=False)

SEQUENCE_LINE=data.readline()
COLUMN_NAME_LINE=data.readline()

COL_NAME_LIST=COLUMN_NAME_LINE.split()
#####################################################################


reweight_info_list=[]
count=0

for reweight_string in reweight_string_list:

	reweight_split=reweight_string.split("-")
	if(len(reweight_split)!=2): error_exit_with_message("len(reweight_split)=(%s)!=2" %(len(reweight_split)))

	reweight_info={}
	weight_name=reweight_split[0]

	try:
		weight_col_index=COL_NAME_LIST.index(weight_name)
	except:
		print "COL_NAME_LIST=", COL_NAME_LIST
		error_exit_with_message("Cannot find column index of weight_name (%s)!" %(weight_name) )


	reweight_factor=reweight_split[1]

	if(reweight_factor.count("/")==0):
		reweight_factor=float(reweight_factor)
	elif(reweight_factor.count("/")==1):
		reweight_factor=float(reweight_factor.split("/")[0])/float(reweight_factor.split("/")[1])
	else:
		error_exit_with_message("reweight_factor.count(\"/\")>1!, reweight_factor=%s" %(reweight_factor))

	if(weight_name=="score"): error_exit_with_message('weight_name=="core""')
	if(weight_name=="description"): error_exit_with_message('weight_name=="description"')
	if(weight_name.count("rms")>0): error_exit_with_message('weight_name.count("rms")>0, weight_name=%s' %(weight_name) )

	reweight_info["col_index"]=weight_col_index
	reweight_info["weight_name"]=weight_name
	reweight_info["reweight_factor"]=reweight_factor

	count+=1

	print "reweight_info #%d: " %(count), reweight_info

	reweight_info_list.append(reweight_info)

############4. Create the new silent_file with the shift_RMSD data add#######

output_silent_file="REWEIGHT"

for reweight_info in reweight_info_list:
	output_silent_file=output_silent_file+"_%s_%0.2f" %(reweight_info["weight_name"], reweight_info["reweight_factor"])

output_silent_file=output_silent_file+"_%s" %(basename(input_silent_file))

print "output_silent_file=%s" %(output_silent_file)



OUTPUT_SILENT_FILE=open(output_silent_file, 'w')


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

 
OUTPUT_SILENT_FILE.write(SEQUENCE_LINE) #Sequence line


#####################################Consistency check################################################
if(COL_NAME_LIST.count("SCORE:")!=1): error_exit_with_message("COL_NAME_LIST.count(\"SCORE:\")!=1")
if(COL_NAME_LIST.index("SCORE:")!=0): error_exit_with_message("COL_NAME_LIST.index(\"SCORE:\")!=0")

if(COL_NAME_LIST.count("description")!=1): error_exit_with_message("COL_NAME_LIST.count(\"description\")!=1")
if(COL_NAME_LIST.index("description")!=(len(COL_NAME_LIST)-1)): 
	error_exit_with_message("COL_NAME_LIST.index(\"description\")!=(len(COL_NAME_LIST)-1)")
######################################################################################################

col_name_line=""

for n in range(len(COL_NAME_LIST)):
	if(n==0):
		col_name_line+="%-*s " %(len(COL_NAME_LIST[n])+2, COL_NAME_LIST[n])
	else:
		col_name_line+=" %*s" %(max(len(COL_NAME_LIST[n])+2,10), COL_NAME_LIST[n])

if(len(col_name_line.split())!=len(COL_NAME_LIST)): 
	error_exit_with_message("len(col_name_line.split())=(%s)!(%s)=len(COL_NAME_LIST)" %(len(col_name_line),len(COL_NAME_LIST)))

OUTPUT_SILENT_FILE.write(col_name_line + '\n')
######################################################################################################

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

		for reweight_info in reweight_info_list:

			old_weight_score=float(cols[reweight_info["col_index"]])
			new_weight_score=reweight_info["reweight_factor"]*old_weight_score

			cols[reweight_info["col_index"]]="%8.3f" %( new_weight_score ) #Dec 10, 2011: USED TO BE "%10.3f" %(total_score) #OK since include extra white_space character below

			total_score=total_score-old_weight_score+new_weight_score

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


