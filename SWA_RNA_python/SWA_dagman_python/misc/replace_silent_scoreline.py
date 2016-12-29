#!/usr/bin/env python

from sys import argv,exit,stdout,stderr
import sys
from os.path import exists
import string

import time

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options

#########################Feb 09, 2012...REBUILD FULL LENGTH REGION AND BULGE  REMOVED some RMSD columns####################################

#replace_silent_scoreline.py -keep_column_names  score     fa_atr     fa_rep    fa_intra_rep    fa_intra_RNA_base_phos_atr    fa_intra_RNA_base_phos_rep    lk_nonpolar    lk_nonpolar_intra_RNA    hack_elec_rna_phos_phos    ch_bond    rna_torsion    rna_sugar_close    fa_stack    hbond_sr_bb_sc    hbond_lr_bb_sc    hbond_sc    hbond_intra    geom_sol_intra_RNA    CI_geom_sol    atom_pair_constraint    angle_constraint    rna_bulge    linear_chainbreak    all_rms       O_rmsd    O_loop_rmsd    O_V_rms    O_V_loop_rms    O_PBP_rmsd  -infile region_3_6_sample.cluster.out


#########################For Comparing CURR AND UPDATED CODE##################################################################################################################
#replace_silent_scoreline.py -keep_column_names score     fa_atr     fa_rep    fa_intra_rep    lk_nonpolar    hack_elec_rna_phos_phos    ch_bond    rna_torsion    rna_sugar_close    fa_stack    hbond_sr_bb_sc    hbond_lr_bb_sc    hbond_sc    CI_geom_sol    rna_bulge    all_rms     O_rmsd    O_loop_rmsd    O_V_rms    O_V_loop_rms    O_PBP_rmsd  -infile cluster_0_region_5_4_sample.out_final_sample


#replace_silent_scoreline.py -keep_column_names score     fa_atr     fa_rep    fa_intra_rep    lk_nonpolar    hack_elec_rna_phos_phos    ch_bond    rna_torsion    rna_sugar_close    fa_stack    hbond_sr_bb_sc    hbond_lr_bb_sc    hbond_sc    CI_geom_sol    atom_pair_constraint    angle_constraint    rna_bulge    linear_chainbreak    all_rms     O_rmsd    O_loop_rmsd    O_V_rms    O_V_loop_rms    O_PBP_rmsd -infile region_0_1_sample.out_before_minimize


############################################################################################################################################################

#replace_silent_scoreline.py -keep_column_names score O_loop_rmsd shift_RMSD shift_score num_shift_data -tag_prestring 1JJ2_   -infile WITH_SHIFT_STATS_rebuild_bulge_clustered_FINAL_filtered_energy.out

#replace_silent_scoreline.py -keep_column_names score O_loop_rmsd shift_RMSD shift_score num_shift_data -tag_prestring RNA_09_ -infile WITH_SHIFT_STATS_rebuild_bulge_clustered_FINAL_filtered_energy.out

#replace_silent_scoreline.py -keep_column_names score O_loop_rmsd shift_RMSD shift_score num_shift_data -tag_prestring SWA_    -infile WITH_SHIFT_STATS_rebuild_bulge_region_FINAL.out 



#replace_silent_scoreline.py -keep_column_names score O_loop_rmsd shift_RMSD shift_score num_shift_data   -infile  Jan_28_WITH_SHIFT_STATS_rebuild_bulge_clustered_FINAL_filtered_energy.out

#replace_silent_scoreline.py -keep_column_names score O_loop_rmsd shift_RMSD shift_score num_shift_data   -infile  March_05_WITH_SHIFT_STATS_rebuild_bulge_ALL_region_FINAL.out

#replace_silent_scoreline.py -keep_column_names score O_loop_rmsd  -infile  rebuild_bulge_region_FINAL.out 

############################################################################################################################################################
#   

#replace_silent_scoreline.py -keep_column_names score O_loop_rmsd shift_RMSD shift_score num_shift_data -tag_prestring 1JJ2_   -infile WITH_SHIFT_STATS_Dec_10_1JJ2_rebuild_bulge_clustered_FINAL_filtered_energy.out 

#replace_silent_scoreline.py -keep_column_names score O_loop_rmsd shift_RMSD shift_score num_shift_data -tag_prestring RNA_09_ -infile WITH_SHIFT_STATS_Dec_10_RNA09_rebuild_bulge_clustered_FINAL_filtered_energy.out

#replace_silent_scoreline.py -keep_column_names score O_loop_rmsd shift_RMSD shift_score num_shift_data -tag_prestring SWA_    -infile WITH_SHIFT_STATS_Dec_10_SWA_rebuild_bulge_region_FINAL.out 

#replace_silent_scoreline.py -keep_column_names score O_loop_rmsd-CRYSTAL_rmsd NEW_O_loop_rmsd-NMR_rmsd shift_RMSD shift_score num_shift_data -infile recal_rmsd_WITH_SHIFT_STATS_combine.out

#replace_silent_scoreline.py -keep_column_names score O_loop_rmsd -tag_prestring RNA_09_ -infile Dec_19_RNA09_rebuild_bulge_clustered_FINAL_filtered_energy.out

#replace_silent_scoreline.py -keep_column_names score O_loop_rmsd -tag_prestring SWA_ -infile Dec_10_SWA_rebuild_bulge_region_FINAL.out  

#replace_silent_scoreline.py -keep_column_names score O_loop_rmsd -tag_prestring 1JJ2_ -infile Dec_19_1JJ2_rebuild_bulge_clustered_FINAL_filtered_energy.out  

############################################################################################################################################################

#replace_silent_scoreline.py -keep_column_names score-Rosetta_energy shift_RMSD-shift_RMSD num_shift_data-num_shift_data -infile select_models_suite_clustering_cluster_200_WITH_SHIFT_STATS_combine.out 

##replace_silent_scoreline.py -keep_column_names score shift_RMSD NEW_Full_L_rmsd -infile recal_full_WITH_SHIFT_STATS_region_FINAL.out

#replace_silent_scoreline.py -keep_column_names score O_loop_rmsd NEW_O_loop_rmsd -tag_prestring 1JJ2_ -infile recal_rmsd_Dec_01_1JJ2_rebuild_bulge_clustered_FINAL_filtered_energy.out

#replace_silent_scoreline.py -keep_column_names score O_loop_rmsd-EXP_1_rmsd NEW_O_loop_rmsd-EXP_2_rmsd shift_RMSD shift_score num_shift_data -infile WITH_SHIFT_STATS_combine.out

#replace_silent_scoreline.py -keep_column_names score O_loop_rmsd NEW_O_loop_rmsd -tag_prestring RNA_09_ -infile recal_rmsd_rebuild_bulge_clustered_FINAL_filtered_energy.out

#replace_silent_scoreline.py -keep_column_names score O_loop_rmsd NEW_O_loop_rmsd -tag_prestring SWA_ -infile recal_rmsd_rebuild_bulge_region_FINAL.out

#replace_silent_scoreline.py -keep_column_names score O_loop_rmsd NEW_O_loop_rmsd -tag_prestring SWA_5K_ -infile recal_rmsd_Dec_01_SWA_5000_rebuild_bulge_region_FINAL.out

#replace_silent_scoreline.py -keep_column_names score O_loop_rmsd shift_RMSD shift_score num_shift_data -tag_prestring SWA_ -infile recal_full_WITH_SHIFT_STATS_rebuild_bulge_region_FINAL.out 


#replace_silent_scoreline.py -keep_column_names score-Rosetta_energy shift_RMSD-shift_RMSD num_shift_data-num_shift_data -infile select_models_cluster_0_WITH_IMINO_PENALTY.out_WITHOUT_IMINO_PENALTY.out

#replace_silent_scoreline.py -keep_column_names score-score O_loop_rmsd O_V_loop_rms -tag_prestring FARFAR_ -infile rebuild_bulge_clustered_FINAL_filtered_energy.out 

#replace_silent_scoreline.py -keep_column_names score-score O_loop_rmsd O_V_loop_rms -tag_prestring SWA_ -infile rebuild_bulge_region_FINAL.out

#replace_silent_scoreline.py -keep_column_names score-Rosetta_energy shift_RMSD-shift_RMSD num_shift_data-num_shift_data -infile select_models_cluster_0_shift_RMSD_abs_score_cut_30_WITH_SHIFT_STATS_combine.out

#-tag_prestring SWA_

#replace_silent_scoreline.py -keep_column_names score-score O_loop_rmsd O_V_loop_rms -tag_prestring FARFAR_ -infile rebuild_bulge_clustered_FINAL_filtered_energy.out 

keep_column_names=parse_options( argv, "keep_column_names", [""])

tag_prestring=parse_options( argv, "tag_prestring", "")

infile=parse_options( argv, "infile", "")

if(keep_column_names==[""]): error_exit_with_message("keep_column_names==[\"\"]")

if(infile==""): error_exit_with_message("infile=\"\"")

if(exists(infile)==False): error_exit_with_message("infile (%s) doesn't exist!")

assert_is_valid_non_empty_silent_file(infile)

outfile=parse_options( argv, "outfile", "")

if(outfile==""):
	outfile="replace_scoreline_%s" %(basename(infile))
else:
	if(exists(outfile)): error_exit_with_message("outfile already exist!")

OUTPUT_SILENT_FILE = open( outfile, "w" )

data=safe_open(infile, mode='r' ,Is_master=False) 

SEQUENCE_LINE = data.readline() # The SEQUENCE: gggcgcagccu line
COLUMN_NAME_LINE   = data.readline() # The column name line
COL_NAME_LIST=COLUMN_NAME_LINE.split()

if(SEQUENCE_LINE[0:9]!='SEQUENCE:'):  		 error_exit_with_message("SEQUENCE_LINE[0:9!='SEQUENCE:' for SEQUENCE_LINE (%s)" %(SEQUENCE_LINE) )
if(COLUMN_NAME_LINE[0:6]!='SCORE:'):  		 error_exit_with_message("COLUMN_NAME_LINE[0:6]!='SCORE:'" )
if(COLUMN_NAME_LINE.find('description') == -1 ):  error_exit_with_message("COLUMN_NAME_LINE.find('description') == -1" )
if(COL_NAME_LIST[0]!='SCORE:'): error_exit_with_message("COL_NAME_LIST[0]!='SCORE:" )
if(COL_NAME_LIST[-1]!='description'): error_exit_with_message("COL_NAME_LIST[-1]!='description'" )

keep_column_names=["SCORE:-SCORE:"] + keep_column_names
keep_column_names.append('description-description')


new_to_old_col_name_map={}

NEW_COL_NAME_LIST=[]


for COL_NAME_pair in keep_column_names:
	
	old_col_name=""
	new_col_name=""

	if(len(COL_NAME_pair.split("-"))==1):
		old_col_name=COL_NAME_pair
		new_col_name=COL_NAME_pair

	elif(len(COL_NAME_pair.split("-"))==2):

		old_col_name=COL_NAME_pair.split("-")[0]
		new_col_name=COL_NAME_pair.split("-")[1]

	else:
		error_exit_with_message("Invalid COL_NAME_pair (%s) " %(COL_NAME_pair))

	if(COL_NAME_LIST.count(old_col_name)!=1): error_exit_with_message("COL_NAME_LIST.count(old_col_name)!=1 for old_col_name=%s" %(old_col_name))

	if(NEW_COL_NAME_LIST.count(new_col_name)!=0): error_exit_with_message("new_col_name (%s) already exist in NEW_COL_NAME_LIST" %(new_col_name))

	NEW_COL_NAME_LIST.append(new_col_name)

	new_to_old_col_name_map[new_col_name]=COL_NAME_LIST.index(old_col_name)

######################################################################################################

OUTPUT_SILENT_FILE.write(SEQUENCE_LINE)

#####################################Consistency check################################################
if(NEW_COL_NAME_LIST.count("SCORE:")!=1): error_exit_with_message("NEW_COL_NAME_LIST.count(\"SCORE:\")!=1")
if(NEW_COL_NAME_LIST.index("SCORE:")!=0): error_exit_with_message("NEW_COL_NAME_LIST.index(\"SCORE:\")!=0")

if(NEW_COL_NAME_LIST.count("description")!=1): error_exit_with_message("NEW_COL_NAME_LIST.count(\"description\")!=1")
if(NEW_COL_NAME_LIST.index("description")!=(len(NEW_COL_NAME_LIST)-1)): 
	error_exit_with_message("NEW_COL_NAME_LIST.index(\"description\")!=(len(NEW_COL_NAME_LIST)-1)")

######################################################################################################
new_col_name_line=""

for col_name in NEW_COL_NAME_LIST:

	if(col_name=="SCORE:"):
		new_col_name_line+="%-*s " %(len(col_name)+2, col_name)
	else:
		new_col_name_line+=" %*s" %(max(len(col_name)+2,10), col_name) 

if(len(new_col_name_line.split())!=len(NEW_COL_NAME_LIST)): 
	error_exit_with_message("len(new_col_name_line.split())=(%s)!(%s)=len(NEW_COL_NAME_LIST)" %(len(new_col_name_line),len(NEW_COL_NAME_LIST)))

OUTPUT_SILENT_FILE.write(new_col_name_line + '\n')

######################################################################################################

for line in data:

	if( (line.find("SCORE:") != -1) != (line[0:6]=="SCORE:") ): error_exit_with_message("(line.find(\"SCORE:\") != -1) != (line[0:6]==\"SCORE:\") for line=%s" %(line))

	if(line[0:6]=="SCORE:"): #SCORE line! 
		
		if(line.find('description') != -1): error_exit_with_message("extra column_name line (%s)" %(line) )

		if(line.find('score') != -1): error_exit_with_message("extra column_name line (%s)" %(line) ) #Feb 06, 2012
		
		score_line=line

		cols=score_line.split()

		if(len(cols)!=len(COL_NAME_LIST)): error_exit_with_message("len(cols)!=len(COL_NAME_LIST)")

		#################################################################################
		new_score_line=""

		for col_name in NEW_COL_NAME_LIST:

			new_col_value=cols[new_to_old_col_name_map[col_name]]

			if(col_name=="SCORE:"):
				new_score_line+="%-*s " %(len(col_name)+2, new_col_value) 

			elif(col_name=="description"):

				new_tag=new_col_value
				if(tag_prestring!=""): new_tag=tag_prestring + new_tag
				new_score_line+=" %*s" %(max(len(col_name)+2,10), new_tag)

			else:
				new_score_line+=" %*s" %(max(len(col_name)+2,10), new_col_value)

		if(len(new_score_line.split())!=len(NEW_COL_NAME_LIST)): 
			error_exit_with_message("len(new_score_line.split())=(%s)!=(%s)=len(NEW_COL_NAME_LIST) for new_score_line=%s" %(len(new_score_line.split()),len(NEW_COL_NAME_LIST), new_score_line))

		OUTPUT_SILENT_FILE.write(new_score_line +'\n')
		#################################################################################

	elif(line[0:7]=="REMARK "): 

		OUTPUT_SILENT_FILE.write(line)

	else:

		description_index = line.rfind(' ')+1 

		tag = line[description_index:-1]	 #-1 to remove the new_line '\n' character

		if(tag!=line.split()[-1]): error_exit_with_message("tag=(%s)!=(%s)=line.split()[-1], line=%s" %(tag, line.split()[-1], line) )

		new_tag=tag
		if(tag_prestring!=""): new_tag=tag_prestring + new_tag

		line = line[:description_index] + new_tag + '\n'

		OUTPUT_SILENT_FILE.write(line)

data.close()

OUTPUT_SILENT_FILE.close()

ensure_no_duplicate_tags(outfile)

