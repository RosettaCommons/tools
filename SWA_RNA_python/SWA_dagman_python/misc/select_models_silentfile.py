#!/usr/bin/env python

from sys import argv,exit,stdout,stderr
import sys
from os.path import exists
import string

import time

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options

######################WARNING, this script will rank models according to order they appear in the silent_file!#####################

###/Dec_19_COMBINE_TRNA_ANTICODON/BLIND_TRNA_A_PLUS_BS_UCC/TRAIL_4_CURR_CODEBASE_CONSISTENCY_CHECK/###
#select_models_silentfile.py   -infile suite_clustering_cluster_200_WITH_SHIFT_STATS_combine.out -model_list S_000000 S_000002 S_000005 S_000009 S_000010 S_000014 S_000017  S_000019  S_000020  S_000021  


###/Dec_19_COMBINE_TRNA_ANTICODON/BLIND_TRNA_BS_UCC/TRAIL_3_3X_CHEM_SHIFT_CURRENT_CODE_BASE/###
#select_models_silentfile.py   -infile suite_clustering_cluster_200_WITH_SHIFT_STATS_combine.out -model_list S_000000 S_000002 S_000003 S_000004 S_000007 S_000012 S_000014  S_000015  S_000019  S_000020

###/Dec_19_COMBINE_TRNA_ANTICODON/BLIND_TRNA_BS_GCC/TRAIL_3_3X_CHEM_SHIFT_SCORE_CURR_BASE_CODE/###
#select_models_silentfile.py   -infile suite_clustering_cluster_200_WITH_SHIFT_STATS_combine.out -model_list S_000000 S_000001 S_000003 S_000004 S_000008 S_000009 S_000020  S_000026  S_000030  S_000031


###/Dec_19_COMBINE_TRNA_ANTICODON/BLIND_TRNA_SE_UCC/TRAIL_4_3X_CHEM_SHIFT_SCORE_CURR_CODE_BASE/###
#select_models_silentfile.py   -infile suite_clustering_cluster_200_WITH_SHIFT_STATS_combine.out -model_list S_000000 S_000001 S_000002 S_000003 S_000005 S_000006 S_000010  S_000013  S_000015  S_000017


####UAAC_tetraloop (WITHOUT_WC)######
#select_models_silentfile.py   -infile suite_clustering_cluster_200_WITH_SHIFT_STATS_combine.out -model_list S_000000 S_000004 S_000007 S_000008 S_000017 S_000019 S_000021  S_000026  S_000027  S_000030 S_000031  S_000032 S_000033  


####UGAC_tetraloop (WITHOUT_WC)######
#select_models_silentfile.py   -infile suite_clustering_cluster_200_WITH_SHIFT_STATS_combine.out -model_list S_000000 S_000001 S_000002 S_000003 S_000004 S_000005 S_000007  S_000008 S_000009  S_000011  S_000016  S_000019  S_000020  S_000024  S_000026



####UUAC_tetraloop (WITH_WC)######
#select_models_silentfile.py   -infile suite_clustering_cluster_200_WITH_SHIFT_STATS_combine.out -model_list S_000000 S_000003 S_000005 S_000007 S_000009 S_000020 S_000022  S_000026 S_000032  S_000038  S_000039  S_000040 


##########UCAC_tetraloop############## 
######(WITH_WC)######
#select_models_silentfile.py   -infile suite_clustering_cluster_200_WITH_SHIFT_STATS_combine.out -model_list S_000000 S_000004 S_000007 S_000011 S_000012 S_000017 S_000021  S_000035 S_000038 -skip_model_rank_list 8 9

######(WITHOUT_WC)######
#select_models_silentfile.py   -infile suite_clustering_cluster_200_WITH_SHIFT_STATS_combine.out -model_list S_000014 S_000017 S_000020 S_000021 S_000023 S_000024 -skip_model_rank_list 1 2 3 4 5 6 7 10 11 -total_score_offset N17.45600

######Alternative M_1 models.
#select_models_silentfile.py   -infile suite_clustering_cluster_200_WITH_SHIFT_STATS_combine.out -model_list S_000001 S_000003 S_000008  -tag_prestring AL_1_

###########Sigel_39mer###########################################

#select_models_silentfile.py   -infile cluster_0_WITH_IMINO_PENALTY.out_WITHOUT_IMINO_PENALTY.out -model_list S_000000 S_000001 S_000004 S_000008 S_000010 S_000016 S_000018  S_000024 S_000006 S_000013 -rank_models True

#select_models_silentfile.py   -infile suite_clustering_cluster_200_REWEIGHT_shift_score_3.00_WITH_SHIFT_STATS_combine.out -model_list S_000000 S_000001 S_000004 S_000008 S_000010 S_000016 S_000018  S_000024 -rank_models False

#select_models_silentfile.py   -infile suite_clustering_cluster_200_REWEIGHT_shift_score_3.00_WITH_SHIFT_STATS_combine.out -model_list S_000006 S_000013 -total_score_offset 10.0 -rank_models False

##################################################################

#CHIMP HAR1F_WITH_LOWER_BP

#select_models_silentfile.py   -infile suite_clustering_cluster_200_shift_RMSD_abs_score_cut_30_WITH_SHIFT_STATS_combine.out -model_list S_000000 S_000001 S_000002 S_000003 S_000005 S_000006 S_000008 S_000009 S_000012 S_000013 S_000014 S_000016 S_000046


#CHIMP HAR1F
#select_models_silentfile.py   -infile cluster_0_shift_RMSD_abs_score_gap_30_WITH_SHIFT_STATS_combine.out -model_list S_0 S_1 S_2 S_7 S_11 S_16 S_36 S_45 S_120

#HUMAN HAR1F
#select_models_silentfile.py   -infile cluster_0_shift_RMSD_abs_score_cut_30_WITH_SHIFT_STATS_combine.out -model_list S_000000 S_000011 S_000027 S_000032 S_000033 S_000035 S_000036 S_000117 S_000120



#S_000000 S_000001 S_000002 S_000007 S_000011 S_000016 S_000036 S_000045 S_000120

#-tag_prestring SWA_

model_list=parse_options( argv, "model_list", [""])

tag_prestring=parse_options( argv, "tag_prestring", "")

rank_models=parse_options( argv, "rank_models", "True")

total_score_offset=parse_options( argv, "total_score_offset", "0.0")

skip_model_rank_list=parse_options( argv, "skip_model_rank_list", [0])

infile=parse_options( argv, "infile", "")

if(infile==""): error_exit_with_message("infile=\"\"")

if(exists(infile)==False): error_exit_with_message("infile (%s) doesn't exist!")

assert_is_valid_non_empty_silent_file(infile)

outfile=parse_options( argv, "outfile", "")

if(outfile==""): outfile="select_models_%s" %(basename(infile))

if(exists(outfile)): 
	print "WARNING outfile (%s) already exist ...removing!" %(outfile)
	submit_subprocess("rm %s" %(outfile))

OUTPUT_SILENT_FILE = open( outfile, "w" )

data=safe_open(infile, mode='r' ,Is_master=False) 

SEQUENCE_LINE = data.readline() # The SEQUENCE: gggcgcagccu line
COLUMN_NAME_LINE   = data.readline() # The column name line
COL_NAME_LIST=COLUMN_NAME_LINE.split()

if(SEQUENCE_LINE[0:9]!='SEQUENCE:'):  		 error_exit_with_message("SEQUENCE_LINE[0:9]!='SEQUENCE:' for SEQUENCE_LINE (%s)" %(SEQUENCE_LINE) )
if(COLUMN_NAME_LINE[0:6]!='SCORE:'):  		 error_exit_with_message("COLUMN_NAME_LINE[0:6]!='SCORE:'" )
if(COLUMN_NAME_LINE.find('description') == -1 ):  error_exit_with_message("COLUMN_NAME_LINE.find('description') == -1" )
if(COL_NAME_LIST[0]!='SCORE:'): error_exit_with_message("COL_NAME_LIST[0]!='SCORE:" )
if(COL_NAME_LIST[-1]!='description'): error_exit_with_message("COL_NAME_LIST[-1]!='description'" )


NEW_COL_NAME_LIST=[]

OUTPUT_SILENT_FILE.write(SEQUENCE_LINE)


#####################################Consistency check################################################
if(COL_NAME_LIST.count("SCORE:")!=1): error_exit_with_message("COL_NAME_LIST.count(\"SCORE:\")!=1")
if(COL_NAME_LIST.index("SCORE:")!=0): error_exit_with_message("COL_NAME_LIST.index(\"SCORE:\")!=0")

if(COL_NAME_LIST.count("description")!=1): error_exit_with_message("COL_NAME_LIST.count(\"description\")!=1")
if(COL_NAME_LIST.index("description")!=(len(COL_NAME_LIST)-1)): 
	error_exit_with_message("COL_NAME_LIST.index(\"description\")!=(len(NEW_COL_NAME_LIST)-1)")

######################################################################################################

col_name_line=""

for col_name in COL_NAME_LIST:
	if(col_name=="SCORE:"):
		col_name_line+="%-*s " %(len(col_name)+2, col_name)
	else:
		col_name_line+=" %*s" %(max(len(col_name)+2,10), col_name)

if(len(col_name_line.split())!=len(COL_NAME_LIST)): 
	error_exit_with_message("len(col_name_line.split())=(%s)!(%s)=len(COL_NAME_LIST)" %(len(col_name_line),len(COL_NAME_LIST)))

OUTPUT_SILENT_FILE.write(col_name_line + '\n')

######################################################################################################

model_count_map={} #Actually could use pigeon hold principle BUT this should be OK

curr_model_rank=0

Is_selected_model=False

for line in data:

	##########################################################################

	if( (line.find("SCORE:") != -1) != (line[0:6]=="SCORE:") ): error_exit_with_message("(line.find(\"SCORE:\") != -1) != (line[0:6]==\"SCORE:\") for line=%s" %(line))
	##########################################################################


	if( (line.split()[0]=="REMARK") and (line.split()[1]=="BINARY_SILENTFILE") and (line.split()[2]=="RNA") ):
		OUTPUT_SILENT_FILE.write(line)
		continue

	##########################################################################
	if(line[0:6]=="SCORE:"): #SCORE line! 
		
		if(line.find('description') != -1): error_exit_with_message("extra column_name line (%s)" %(line) )

		score_line=line

		cols=score_line.split()

		if(len(cols)!=len(COL_NAME_LIST)): error_exit_with_message("len(cols)!=len(COL_NAME_LIST)")

		old_tag=cols[-1]

		if(model_count_map.has_key(old_tag)==True): error_exit_with_message("old_tag (%s) already exist in model_count_map!" %(old_tag))

		if old_tag in model_list:
			model_count_map[old_tag]=True
			curr_model_rank+=1
			Is_selected_model=True
		else:
			Is_selected_model=False

		if(len(skip_model_rank_list)!=0):
			while(curr_model_rank in skip_model_rank_list):  
				curr_model_rank+=1

		if(rank_models):
			
			if(old_tag.count("_")==0): error_exit_with_message("old_tag.count(\"_\")==0")

			start_rank_index= old_tag.rfind('_')+1 

			new_tag = "%s%s" %(old_tag[:start_rank_index] , curr_model_rank)

		else:
			new_tag=old_tag

		if(tag_prestring!=""): new_tag=tag_prestring + new_tag



	############################################################################

	if(Is_selected_model==False): continue

	if(line[0:6]=="SCORE:"): #SCORE line! 

		if(line.find('description') != -1): error_exit_with_message("extra column_name line (%s)" %(line) )
		
		score_line=line

		cols=score_line.split()

		if(len(cols)!=len(COL_NAME_LIST)): error_exit_with_message("len(cols)!=len(COL_NAME_LIST)")

		####################################################################################################
		new_score_line=""

		for col_num in range(len(COL_NAME_LIST)):

			col_value=cols[col_num]
			col_name=COL_NAME_LIST[col_num]

			if(col_name=="SCORE:"):
				new_score_line+="%-*s " %(len(col_name)+2, col_value)

			elif(col_name=="score"):
				if(total_score_offset=="0.0"):
					score_value=col_value
				else:
					total_scrore_offset_float=convert_string_to_float(total_score_offset)

					score_value="%8.3f" %(float(col_value)+total_scrore_offset_float)

				new_score_line+=" %*s" %(max(len(col_name)+2,10), score_value)

			elif(col_name=="description"):

				new_score_line+=" %*s" %(max(len(col_name)+2,10), new_tag)

			else:
				new_score_line+=" %*s" %(max(len(col_name)+2,10), col_value)

		if(len(score_line.split())!=len(COL_NAME_LIST)): 
			error_exit_with_message("len(new_score_line.split())=(%s)!=(%s)=len(COL_NAME_LIST) for new_score_line=%s" %(len(new_score_line.split()),len(COL_NAME_LIST), new_score_line))
		######################################################################################################

		OUTPUT_SILENT_FILE.write(new_score_line +'\n')

	elif(line[0:7]=="REMARK "): 

		OUTPUT_SILENT_FILE.write(line)

	else:

		description_index = line.rfind(' ')+1 

		tag = line[description_index:-1]	 #-1 to remove the new_line '\n' character

		if(tag!=line.split()[-1]): error_exit_with_message("tag=(%s)!=(%s)=line.split()[-1], line=%s" %(tag, line.split()[-1], line) )

		line = line[:description_index] + new_tag + '\n'

		OUTPUT_SILENT_FILE.write(line)

data.close()

OUTPUT_SILENT_FILE.close()

ensure_no_duplicate_tags(outfile)

########################################################################################################

if((curr_model_rank-len(skip_model_rank_list))!=len(model_list)): error_exit_with_message("(curr_model_rank-len(skip_model_rank_list))=(%s)!=(%s)=len(model_list)" %((curr_model_rank-len(skip_model_rank_list)),len(model_list)) )

for model_name in model_list:

	if(model_count_map.has_key(model_name)==False): error_exit_with_message("model_name (%s) doesn't exist in model_count_map!" %(model_name))

	if(model_count_map[model_name]!=True): error_exit_with_message("model_count_map[model_name]!=True | model_name=%s" %(model_name))




