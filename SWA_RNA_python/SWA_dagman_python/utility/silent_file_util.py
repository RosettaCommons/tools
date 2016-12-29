#!/usr/bin/env python

from error_util import *
#######################################################
'''This function have the same name will the one in DAG_rebuild_bulge_util.py ...should just integrate them together! them!
def parse_annotated_sequence(ANNOTATED_SEQUENCE):

	annotated_sequence_data={}

	seq_num=0
	seq_list=[]
	virtual_res_list=[]
	cutpoint_lower_list=[]
	cutpoint_upper_list=[]
	virtual_ribose_list=[]
	n=0

	#ANNOTATED_SEQUENCE: g[RGU_p:Virtual_Phosphate]ccga[RAD_p:Virtual_RNA_Residue]a[RAD_p:rna_cutpoint_lower_p:Virtual_Phosphate_p:Virtual_RNA_Residue_Upper]a[RAD_p:rna_cutpoint_upper]ggc S_7
	while(n<len(ANNOTATED_SEQUENCE)):

		char=ANNOTATED_SEQUENCE[n]

		if(char=='['): 
			end_bracket=ANNOTATED_SEQUENCE.find(']', n)	
	
			Is_virtual_res=False

      #Need to be careful to distinguish Virtual_RNA_Residue and Virtual_RNA_Residue_Upper!
			if( ANNOTATED_SEQUENCE.count('Virtual_RNA_Residue:', n, end_bracket)>0  ): Is_virtual_res=True

			if( ANNOTATED_SEQUENCE.count('Virtual_RNA_Residue]', n, end_bracket+1)>0 ): Is_virtual_res=True

			if( ANNOTATED_SEQUENCE.count('Virtual_RNA_Residue_p:', n, end_bracket)>0  ): Is_virtual_res=True

			if( ANNOTATED_SEQUENCE.count('Virtual_RNA_Residue_p]', n, end_bracket+1)>0 ): Is_virtual_res=True

			Is_cutpoint_lower=(ANNOTATED_SEQUENCE.count('rna_cutpoint_lower', n, end_bracket+1)>0 ) #something weird rna_cutpoint_lower_p (sometime have _p, sometime don't)
			Is_cutpoint_upper=(ANNOTATED_SEQUENCE.count('rna_cutpoint_upper', n, end_bracket+1)>0 )
			Is_virtual_ribose=(ANNOTATED_SEQUENCE.count('Virtual_Ribose', n, end_bracket+1)>0 )

			if(Is_virtual_res): virtual_res_list.append(seq_num)
			if(Is_cutpoint_lower): cutpoint_lower_list.append(seq_num)
			if(Is_cutpoint_upper): cutpoint_upper_list.append(seq_num)
			if(Is_virtual_ribose): virtual_ribose_list.append(seq_num)

			n=end_bracket+1
		else:
			seq_list.append(char )
			seq_num=seq_num+1
			n=n+1
			
	annotated_sequence_data["virtual_res_list"]=virtual_res_list
	annotated_sequence_data["cutpoint_lower_list"]=cutpoint_lower_list
	annotated_sequence_data["cutpoint_upper_list"]=cutpoint_upper_list
	annotated_sequence_data["virtual_ribose_list"]=virtual_ribose_list

	#annotated_sequence_data["seq_list"]=seq_list
	annotated_sequence_data["total_residue"]=seq_num

	#if( len(annotated_sequence_data["seq_list"]) != annotated_sequence_data["total_residue"] ):
	#	error_exit_with_message('len(annotated_sequence_data["seq_list"]) != annotated_sequence_data"total_residue"]')

	return (annotated_sequence_data)
'''

#########################################Copy from SWA_extract_pdb.py#####################################################################################
def get_description_column(silent_file, verbose=True):

	if(exists(silent_file)==False): error_exit_with_message("silent_file (%s) doesn't exist!" (silent_file) )

	description_line=""

	found_description_line=False

	for line in open(silent_file):
	
		if(line.find('description') != -1): #found the description line

			if( (line.find("SCORE") != -1) != (line[0:6]=="SCORE:") ): error_exit_with_message("(line.find(\"SCORE\") != -1) != (line[0:6]==\"SCORE:\")")

			if(line[0:6]!="SCORE:"): error_exit_with_message("line.find('desc') != -1 but line[0:6]!=/SCORE:/")

			found_description_line=True  
			description_line=line
			break;

	if( found_description_line == False ): error_exit_with_message( "found_description_line == False" )

	col_name_list=description_line.split()

	try:
		desc_col_index=col_name_list.index('description')
	except:
		error_exit_with_message("Cannot find score_col_index with corresponding score_colname= %s " %(scorecol_name) )


	if(verbose): print "description is located at column_num: %d" %(desc_col_index+1)
	
	return desc_col_index

