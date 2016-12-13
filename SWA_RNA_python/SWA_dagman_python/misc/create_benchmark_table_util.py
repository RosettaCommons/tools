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
from os import popen ##don't COPY THIS!
######################################################################

from SWA_dagman_python.misc.SWA_cat_outfiles import concatenate_outfiles

from SWA_dagman_python.utility.extract_MC_annotate_data_functions import extract_MC_annotate_data_func
from SWA_dagman_python.utility.extract_FR3D_data_functions import extract_FR3D_data_func
from SWA_dagman_python.utility.RNA_BP_and_BS_util import *
from SWA_dagman_python.utility.DAGMAN_util import parse_chemical_shift_args
from SWA_dagman_python.utility.SWA_util import *

from SWA_dagman_python.parser.SWA_parse_options import parse_options, get_option_name_args_safe
from SWA_dagman_python.parser.SWA_parse_internal_arguments import *
from SWA_dagman_python.parser.SWA_parse_benchmark import parse_benchmark_table_data_file
######################################################################

###Global variables######
print "----------------------------------GLOBAL VARIABLES:----------------------------------"


RMSD_COL_NAME_PREFIX="RMSD_V1" #THIS IS NEW_Full_L_rmsd, (RMSD_V2 is NEW_NAT_rmsd)

RMSD_TO_NAT_COL_NAME="%s_TO_NAT" %(RMSD_COL_NAME_PREFIX)  

NO_GRAPHIC=False
PERFORM_VDW_REP_SCREEN=True

TOP_ENERGY_CLUSTERS_FOLDER="TOP_ENERGY_CLUSTERS/"

OUTPUT_TOP_ENERGY_CLUSTERS_MODE=False
CHEM_SHIFT_DATA_MODE=True

CONFIDENCE_GAP_RMSD=2.0

NATIVE_RMSD_CUTOFF= 1.5
NUM_ENERGY_CLUSTERS= 5  
CLUSTER_RMSD=1.5 
SUITE_CLUSTER_RMSD= 2.5 
DOUBLE_LINE_SPACING=True
DISTINGUISH_CANONICAL_BP=True

DOUBLE_COUNT_BP_AND_BS=True #WARNING, still always double_count for Rosetta BP!
INCLUDE_BIFURCARCATION_BP=True

IGNORE_FARFAR_NO_AUTO_BULGE_PARENT_TAG=False
IGNORE_UNMATCHED_VIRTUAL_RES=False

#############Nov 05, 2011 Chemical_shift stuff#############
OUTPUT_CHEM_SHIFT= False

###########################################################

if(CHEM_SHIFT_DATA_MODE):
	col_type_list=['best_ENERGY', 'best_shift_RMSD']	
	subcol_type_list=["| energy", "rms_NAT", "rms_Bshift", "rms_Benergy", "shift_RMSD", "shift_gapC", "Ros_E_gapC", "E_gapC"]


elif(OUTPUT_TOP_ENERGY_CLUSTERS_MODE):
	col_type_list=[]	
	for n in range(NUM_ENERGY_CLUSTERS): col_type_list.append("cluster_%d" %(n+1) )
	subcol_type_list=["energy", "rmsd"]
	if(OUTPUT_CHEM_SHIFT): subcol_type_list.extend(['shift_RMSD', 'shift_score'])
else:
	#col_type_list=["best_score", "lowest_rmsd", "best_RMSD_cluster", "best_ENERGY_cluster" ]
	#Also note implementation below REQUIRES THAT "best_RMSD_cluster" COME RIGHT BEFORE "best_ENERGY_cluster" (BUT OK it either or both is missing)! 
	#col_type_list=["lowest_rmsd", "best_ENERGY_cluster" ]
	col_type_list=["best_score", "lowest_rmsd"]
	subcol_type_list=["energy", "rmsd"]
	if(OUTPUT_CHEM_SHIFT): subcol_type_list.extend(['shift_RMSD', 'shift_score'])

DO_output_Rosetta_BP_stats=False
DO_output_MC_stats=False
DO_output_FR3D_stats=False



###########The order the subcol_type are add does matter!#################

if(DO_output_Rosetta_BP_stats): subcol_type_list.append("Ros_Pair")
if(DO_output_MC_stats):         subcol_type_list.extend(["MC_Pair","MC_Stack"]) 
if(DO_output_FR3D_stats): 			subcol_type_list.extend(["FR3D_Pair","FR3D_Stack"]) 

INDIV_COL_WIDTH=12
SUB_COLUMNS_WIDTH=len(subcol_type_list)*INDIV_COL_WIDTH
COLUMNS_WIDTH=len(col_type_list)*SUB_COLUMNS_WIDTH

print "--------------------------------------------------------------------"
print "CONFIDENCE_GAP_RMSD= %s " %(CONFIDENCE_GAP_RMSD)
print "CLUSTER_RMSD= %s " %(CLUSTER_RMSD) 
print "SUITE_CLUSTER_RMSD= %s " %(SUITE_CLUSTER_RMSD)
print "NATIVE_RMSD_CUTOFF= %s " %(NATIVE_RMSD_CUTOFF) 
print "--------------------------------------------------------------------"
print "OUTPUT_TOP_ENERGY_CLUSTERS_MODE=%s" %(OUTPUT_TOP_ENERGY_CLUSTERS_MODE)
print "NUM_ENERGY_CLUSTERS= %d " %(NUM_ENERGY_CLUSTERS) 
print "--------------------------------------------------------------------"
print "DOUBLE_COUNT_BP_AND_BS=%s" %(DOUBLE_COUNT_BP_AND_BS)
print "INCLUDE_BIFURCARCATION_BP=%s" %(INCLUDE_BIFURCARCATION_BP)
print "--------------------------------------------------------------------"
print "IGNORE_FARFAR_NO_AUTO_BULGE_PARENT_TAG=%s" %(IGNORE_FARFAR_NO_AUTO_BULGE_PARENT_TAG)
print "IGNORE_UNMATCHED_VIRTUAL_RES=%s" %(IGNORE_UNMATCHED_VIRTUAL_RES)
print "--------------------------------------------------------------------"
print "OUTPUT_CHEM_SHIFT=%s " %(OUTPUT_CHEM_SHIFT)
print "--------------------------------------------------------------------"
print "RMSD_COL_NAME_PREFIX= %s " %(RMSD_COL_NAME_PREFIX)
print "--------------------------------------------------------------------"


print "----------------------------------GLOBAL VARIABLES:----------------------------------"

######################################################################################################

def print_column_names(TABLE, algorithm_list):


	####1ST LINE#####
	TABLE.write( '%-20s' %('Motif_Name') )

	TABLE.write( '%-*s' %(int(COLUMNS_WIDTH/2)+9, '') )

	for algorithm in algorithm_list:		
		TABLE.write( '%-*s' %(COLUMNS_WIDTH, algorithm) ) #2*12*3 

	TABLE.write('\n')

	####2ND LINE#####
	TABLE.write( '%-*s' %(20+int(SUB_COLUMNS_WIDTH/2)+8, '') )

	#TABLE.write( '%-14s' %('') )

	for algorithm in algorithm_list:
		for col_type in col_type_list:	
			TABLE.write( '%-*s' %(SUB_COLUMNS_WIDTH, col_type) )
	TABLE.write('\n')

	####3RD LINE#####
	TABLE.write( '%-*s' %(20+int(INDIV_COL_WIDTH/2)+8, '') ) 
	#TABLE.write( '%-34s' %('') ) #20 +14?
	for algorithm in algorithm_list:
		for n in range(len(col_type_list)):	
			for ii in range(len(subcol_type_list)):	
				TABLE.write( '%-*s' %(INDIV_COL_WIDTH, subcol_type_list[ii]) )
	TABLE.write('\n')


##############################################################################################################################



#####################################################################################################

def get_col_index_value(col_name, data_list, col_name_list):

	try:
		col_index=col_name_list.index(col_name)
	except:
		print "col_name_list=", col_name_list
		error_exit_with_message('Cannot find  score_col= "%s" ' %(col_name))

	return col_index

######################################################################################################

def get_col_score_value(col_name, data_list, col_name_list):

	if(len(data_list)==0): return 99.99
	if(len(col_name_list)==0): return 99.99

	col_index=get_col_index_value(col_name, data_list, col_name_list)

	return float(data_list[col_index])



######################################################################################################
def output_FR3D_stats_func(LOCAL_TABLE, output_foldername, model_pdb, native_pdb, sample_segment, virtual_res_list):

	BASE_folder=os.getcwd()

	FR3D_folder="%s/extract_MC_data/" %(output_foldername)
	submit_subprocess("mkdir -p %s/" %(FR3D_folder) )
	submit_subprocess("mkdir -p %s/LOG/" %(FR3D_folder) )

	os.chdir( FR3D_folder )

	(model_base_pair_list, model_base_stack_list)  = \
		extract_FR3D_data_func( model_pdb , sample_segment, double_count=DOUBLE_COUNT_BP_AND_BS, allow_near_match=False, Verbose=False)

	(native_base_pair_list, native_base_stack_list)= \
		extract_FR3D_data_func( native_pdb, sample_segment, double_count=DOUBLE_COUNT_BP_AND_BS, allow_near_match=False, Verbose=False)

	os.chdir( BASE_folder)

	(num_native_WC, num_native_NWC, recovered_WC, recovered_NWC)=get_base_pair_recovery_stat( native_base_pair_list, model_base_pair_list , virtual_res_list, Verbose=False)

	(num_native_BS, recovered_BS)=get_base_stack_recovery_stat( native_base_stack_list, model_base_stack_list , virtual_res_list, Verbose=False)

	if(DISTINGUISH_CANONICAL_BP):
		LOCAL_TABLE.write( '%-12s' %( '%s/%s|%s/%s' %(recovered_WC, num_native_WC, recovered_NWC, num_native_NWC ) ) )
	else:
		LOCAL_TABLE.write( '%-12s' %( '%s/%s' %(recovered_WC+recovered_NWC, num_native_WC+num_native_NWC ) ) )

	LOCAL_TABLE.write( '%-12s' %( '%s/%s' %(recovered_BS,num_native_BS  ) ) )


######################################################################################################
def output_MC_stats_func(LOCAL_TABLE, output_foldername, model_pdb, native_pdb, sample_segment, virtual_res_list):

	BASE_folder=os.getcwd()

	MC_folder="%s/extract_MC_data/" %(output_foldername)
	submit_subprocess("mkdir -p %s/" %(MC_folder) )
	submit_subprocess("mkdir -p %s/LOG/" %(MC_folder) )

	os.chdir( MC_folder )

	(model_base_pair_list, model_base_stack_list)  = \
		extract_MC_annotate_data_func( model_pdb , sample_segment, double_count=DOUBLE_COUNT_BP_AND_BS, include_bifurcation_BP=INCLUDE_BIFURCARCATION_BP, Verbose=False)

	(native_base_pair_list, native_base_stack_list)= \
		extract_MC_annotate_data_func( native_pdb, sample_segment, double_count=DOUBLE_COUNT_BP_AND_BS, include_bifurcation_BP=INCLUDE_BIFURCARCATION_BP, Verbose=False)

	os.chdir( BASE_folder)

	(num_native_WC, num_native_NWC, recovered_WC, recovered_NWC)=get_base_pair_recovery_stat( native_base_pair_list, model_base_pair_list , virtual_res_list, Verbose=False)

	(num_native_BS, recovered_BS)=get_base_stack_recovery_stat( native_base_stack_list, model_base_stack_list , virtual_res_list, Verbose=False)

	if(DISTINGUISH_CANONICAL_BP):
		LOCAL_TABLE.write( '%-12s' %( '%s/%s|%s/%s' %(recovered_WC, num_native_WC, recovered_NWC, num_native_NWC ) ) )
	else:
		LOCAL_TABLE.write( '%-12s' %( '%s/%s' %(recovered_WC+recovered_NWC, num_native_WC+num_native_NWC ) ) )

	LOCAL_TABLE.write( '%-12s' %( '%s/%s' %(recovered_BS,num_native_BS  ) ) )

######################################################################################################
def output_third_party_stats( LOCAL_TABLE, spec_job, data_list, col_name_list, silent_file, silent_file_desc, recreate_file):

	native_pdb=spec_job["native_pdb"]

	native_pdb=os.path.abspath(native_pdb)

	common_args_file=spec_job["common_args"]

	extract_pdb_folder="%s/EXTRACTED_PDB/" %(TOP_ENERGY_CLUSTERS_FOLDER)

	submit_subprocess("mkdir -p %s" %(extract_pdb_folder) )

	if(silent_file==""): error_exit_with_message("global_silent_file==\"\"" %(silent_file) )
	if(basename(silent_file)==""): error_exit_with_message("basename(silent_file)==\"\" for silent_file (%s)" %(silent_file))

	if(len(data_list)==0): # Not valid data. EARLY RETURN!
		if(DO_output_MC_stats):   LOCAL_TABLE.write( '%-12s%-12s' %( '%s/%s|%s/%s' %("E", "E", "E", "E"), '%s/%s' %("E", "E") ) ) #Output BP and BS STATS
		if(DO_output_FR3D_stats): LOCAL_TABLE.write( '%-12s%-12s' %( '%s/%s|%s/%s' %("E", "E", "E", "E"), '%s/%s' %("E", "E") ) ) #Output BP and BS STATS

		return
		

	if(len(col_name_list)==0): error_exit_with_message("len(col_name_list)==0")

	tag=data_list[get_col_index_value('description', data_list, col_name_list)]

	output_foldername="%s/%s_%s/" %(extract_pdb_folder, silent_file_desc, basename( silent_file ) )

	output_pdb=os.path.abspath(output_foldername + tag + ".pdb")

	############################################################################################################
	if(recreate_file):
		if(exists(output_pdb)):
			print "Warning...recreate_file=True...removing output_pdb (%s) " %(output_pdb) 
			submit_subprocess("rm %s " %(output_pdb) )

	if(exists(output_pdb)==False):
		extract_command="SWA_extract_pdb.py -silent_file %s -output_foldername %s -tag %s " %(os.path.abspath(silent_file), output_foldername, tag)
		extract_command+="-no_graphic %s " %(NO_GRAPHIC)
		submit_subprocess( extract_command )

	if(exists(output_pdb)==False): error_exit_with_message("output_pdb doesn't exist! output_pdb = %s" %( output_pdb ) )

	############################################################################################################

	common_args_list=safe_readlines( common_args_file )[0].split()

	
	sample_res_list =  parse_options(common_args_list, "rmsd_res",            [0], Verbose=False) #Should change to global_sample_res_list!

	virtual_res_list = parse_options(common_args_list, "native_virtual_res", [0], Verbose=False) #Should change to global_sample_res_list!

	#print "native_pdb=%s | sample_res_list=%s | virtual_res_list=%s " %(basename(native_pdb), sample_res_list, virtual_res_list)

	start_res=sample_res_list[0]
	end_res  =sample_res_list[-1] 
	sample_segment="%s-%s" %( start_res, end_res ) 

	#############Consistency check!###################
	for n in range(len(sample_res_list)):
		if(sample_res_list[n]!=start_res+n): error_exit_with_message("sample_res_list[n]!=start_res+n, sample_res_list[n]=%d, start_res=%d, n=%d" %( sample_res_list[n], start_res, n ) )		
	##################################################

	if(DO_output_MC_stats):   output_MC_stats_func(  LOCAL_TABLE, output_foldername, output_pdb, native_pdb, sample_segment, virtual_res_list)

	if(DO_output_FR3D_stats): output_FR3D_stats_func(LOCAL_TABLE, output_foldername, output_pdb, native_pdb, sample_segment, virtual_res_list)



######################################################################################################


def output_Rosetta_BP_stats( LOCAL_TABLE, local_data_list, local_col_name_list):
	
	#NAT_WC    NAT_NWC     REC_WC    REC_NWC

	num_native_WC =  int(get_col_score_value("NAT_WC", local_data_list, local_col_name_list))
	num_native_NWC=  int(get_col_score_value("NAT_NWC", local_data_list, local_col_name_list))
	recovered_WC  = 	int(get_col_score_value("REC_WC", local_data_list, local_col_name_list))
	recovered_NWC =	int(get_col_score_value("REC_NWC", local_data_list, local_col_name_list))

	if(num_native_WC==99):  num_native_WC = 'E'
	if(num_native_NWC==99): num_native_NWC= 'E'
	if(recovered_WC==99):   recovered_WC  = 'E'
	if(recovered_NWC==99):  recovered_NWC = 'E'

	if(DISTINGUISH_CANONICAL_BP):
		LOCAL_TABLE.write( '%-12s' %( '%s/%s|%s/%s' %(recovered_WC, num_native_WC, recovered_NWC, num_native_NWC ) ) )
	else:
		LOCAL_TABLE.write( '%-12s' %( '%s/%s' %(recovered_WC+recovered_NWC, num_native_WC+num_native_NWC ) ) )


##############################################################################################################################
def recalculate_silent_file(input_silent_file_list, native_pdb, common_args_file, recal_rmsd_silent_file): #Checked on March 10, 2012

	print
	print "recalculate_silent_file FOR: input_silent_file_list=(%s) | native_pdb=(%s)" %(list_to_string(input_silent_file_list), basename(native_pdb))


	for input_silent_file in input_silent_file_list:
		assert_is_valid_non_empty_silent_file(input_silent_file)

	if(exists(recal_rmsd_silent_file)): error_exit_with_message("recal_rmsd_silent_file (%s) already exist!" %(recal_rmsd_silent_file))

	recal_rmsd_command="SWA_cluster.py  -num_pose_kept 9999999 -skip_clustering True -recreate_silent_struct True -extract_pdb False "

	recal_rmsd_command+="-write_score_only False "

	recal_rmsd_command+="-silent_file %s "%(list_to_string(input_silent_file_list))
	recal_rmsd_command+="-native_pdb %s "%(native_pdb)
	recal_rmsd_command+="-common_args %s "%(common_args_file)
	recal_rmsd_command+="-output_filename %s " %(recal_rmsd_silent_file)
	recal_rmsd_command+="-perform_VDW_rep_screen %s " %(PERFORM_VDW_REP_SCREEN)
	recal_rmsd_command+="-no_graphic %s " %(NO_GRAPHIC)
	recal_rmsd_command+="-clusterer_rename_tags false " 
	recal_rmsd_command+=" > LOG_recalculate_silent_file_%s.txt" %(basename(native_pdb))	


	submit_subprocess(recal_rmsd_command)

	assert_is_valid_non_empty_silent_file(recal_rmsd_silent_file)



##################################################################################################################################

def create_cluster_silent_file(native_pdb, common_args_file, silent_file_list, recreate_file):


	for silent_file in silent_file_list:
		assert_is_valid_non_empty_silent_file(silent_file)

	cluster_silent_file="%s/top_energy_clusters.out" %(TOP_ENERGY_CLUSTERS_FOLDER)

	if(recreate_file and exists(cluster_silent_file)==True):
		print "recreate_file==True, cluster_silent_file (%s) already exist... removing..." %(os.path.abspath(cluster_silent_file))
		submit_subprocess('rm %s' %(cluster_silent_file))		

	if(exists(cluster_silent_file)==False):
		print
		print "create_cluster_silent_file FOR: cluster_silent_file=(%s)" %(cluster_silent_file)


		submit_subprocess( "mkdir -p %s " %(TOP_ENERGY_CLUSTERS_FOLDER ) )

		#-write_score_only True only work in recalculate RMSD mode...

		cluster_top_energy_command="SWA_cluster.py  -num_pose_kept 100 -distinguish_pucker false -extract_pdb False " 
		cluster_top_energy_command+="-cluster_rmsd %s " %(CLUSTER_RMSD)
		cluster_top_energy_command+="-suite_cluster_rmsd %s " %(SUITE_CLUSTER_RMSD)
		cluster_top_energy_command+="-silent_file %s "%(list_to_string(silent_file_list))
		cluster_top_energy_command+="-native_pdb %s "%(native_pdb)
		cluster_top_energy_command+="-common_args %s "%(common_args_file)
		cluster_top_energy_command+="-output_filename %s " %(cluster_silent_file)	
		cluster_top_energy_command+="-no_graphic %s " %(NO_GRAPHIC)
		cluster_top_energy_command+="-full_length_loop_rmsd_clustering True "
		cluster_top_energy_command+=" > LOG_create_cluster_silent_file.txt"	

		if(IGNORE_FARFAR_NO_AUTO_BULGE_PARENT_TAG):
			cluster_top_energy_command+="-ignore_FARFAR_no_auto_bulge_parent_tag True "

		if(IGNORE_UNMATCHED_VIRTUAL_RES):
			cluster_top_energy_command+="-ignore_unmatched_virtual_res True "

	########################################################################################################
		submit_subprocess(cluster_top_energy_command)
	else:
		print
		print "cluster_silent_file (%s) already exist!" %(cluster_silent_file)

	assert_is_valid_non_empty_silent_file(cluster_silent_file)

	return cluster_silent_file

##################################################################################################################################

def output_top_energy_clusters_stats(TABLE, SCORE_cluster_silent_file, cluster_silent_file, spec_job,  recreate_file):

	if(exists(SCORE_cluster_silent_file)==False): error_exit_with_message("SCORE_cluster_silent_file (%s) doesn't exists!" %(SCORE_cluster_silent_file) )

	common_args_file=spec_job["common_args"]

	lines=open(SCORE_cluster_silent_file).readlines()

	col_name_list=lines[0].split()	

	num_structs=0

	for line_ID in range(1,len(lines)):

		if(line_ID>NUM_ENERGY_CLUSTERS): break

		num_structs+=1

		data_list=lines[line_ID].split()			

		rmsd=float(get_col_score_value(RMSD_TO_NAT_COL_NAME, data_list, col_name_list))
		score=float(get_col_score_value('score', data_list, col_name_list))

		TABLE.write( '%-12s' %(score) )
		TABLE.write( '%-12s' %(rmsd)  )

		if(OUTPUT_CHEM_SHIFT): 
			TABLE.write( '%-12s' %(float(get_col_score_value('R_shift_RMSD', data_list, col_name_list))) )
			TABLE.write( '%-12s' %(float(get_col_score_value('rna_chem_shift', data_list, col_name_list))) )


		if(DO_output_Rosetta_BP_stats): output_Rosetta_BP_stats( TABLE, data_list, col_name_list)

		if(DO_output_MC_stats or DO_output_FR3D_stats): 
			output_third_party_stats( TABLE, spec_job, data_list, col_name_list, cluster_silent_file, "ENERGY_CLUSTER_#%d" %(line_ID), recreate_file)


	if(num_structs<NUM_ENERGY_CLUSTERS):
		error_exit_with_message("num_structs<NUM_ENERGY_CLUSTERS!, num_structs=%s, NUM_ENERGY_CLUSTERS=%d" %(num_structs, NUM_ENERGY_CLUSTERS ) )



##################################################################################################################################

def output_best_energy_and_best_RMSD_cluster_stats(TABLE, SCORE_cluster_silent_file, cluster_silent_file, spec_job, col_type_list, recreate_file):

	if(exists(SCORE_cluster_silent_file)==False): error_exit_with_message("SCORE_cluster_silent_file (%s) doesn't exists!" %(SCORE_cluster_silent_file) )

	common_args_file=spec_job["common_args"]

	best_RMSD_cluster_rmsd=99.0
	best_RMSD_cluster_score=999.0

	best_ENERGY_cluster_rmsd=99.0
	best_ENERGY_cluster_score=999.0

	best_RMSD_cluster_rank=0
	best_ENERGY_cluster_rank=0

	best_RMSD_data_list=[]
	best_RMSD_col_name_list=[]

	best_ENERGY_data_list=[]
	best_ENERGY_col_name_list=[]


	lines=open(SCORE_cluster_silent_file).readlines()

	col_name_list=lines[0].split()	

	num_structs=0

	for line_ID in range(1,len(lines)):

		num_structs+=1
		data_list=lines[line_ID].split()			

		rmsd=float(get_col_score_value(RMSD_TO_NAT_COL_NAME, data_list, col_name_list))
		score=float(get_col_score_value('score', data_list, col_name_list))

		if((rmsd<best_ENERGY_cluster_rmsd) and line_ID<=NUM_ENERGY_CLUSTERS):
			best_ENERGY_cluster_rank=line_ID
			best_ENERGY_cluster_rmsd=rmsd
			best_ENERGY_cluster_score=score
			best_ENERGY_data_list=copy.deepcopy(data_list)
			best_ENERGY_col_name_list=copy.deepcopy(col_name_list)

		if((score<best_RMSD_cluster_score) and rmsd<=NATIVE_RMSD_CUTOFF):
			best_RMSD_cluster_rank=line_ID
			best_RMSD_cluster_rmsd=rmsd
			best_RMSD_cluster_score=score
			best_RMSD_data_list=copy.deepcopy(data_list)
			best_RMSD_col_name_list=copy.deepcopy(col_name_list)

	not_enough_cluster=False
	if(num_structs<NUM_ENERGY_CLUSTERS): not_enough_cluster=True


	if(col_type_list.count("best_RMSD_cluster")==1):
		TABLE.write( '%-12s' %(best_RMSD_cluster_score) )
		TABLE.write( '%-12s' %('%s(%d)' %(best_RMSD_cluster_rmsd, best_RMSD_cluster_rank) ) )

		if(OUTPUT_CHEM_SHIFT): 
			TABLE.write( '%-12s' %(float(get_col_score_value('R_shift_RMSD', best_RMSD_data_list, best_RMSD_col_name_list))) )
			TABLE.write( '%-12s' %(float(get_col_score_value('rna_chem_shift', best_RMSD_data_list, best_RMSD_col_name_list))) )


		if(DO_output_Rosetta_BP_stats): output_Rosetta_BP_stats( TABLE, best_RMSD_data_list, best_RMSD_col_name_list)

		if(DO_output_MC_stats or DO_output_FR3D_stats): 
			output_third_party_stats( TABLE, spec_job, best_RMSD_data_list, best_RMSD_col_name_list, cluster_silent_file, "BEST_RMSD_CLUSTER",  recreate_file)

	if(col_type_list.count("best_ENERGY_cluster")==1):
		TABLE.write( '%-12s' %(best_ENERGY_cluster_score) )
		if(not_enough_cluster):
			TABLE.write( '%-12s' %('%s(!%d)' %(best_ENERGY_cluster_rmsd, best_ENERGY_cluster_rank) ) )
		else:
			TABLE.write( '%-12s' %('%s(%d)' %(best_ENERGY_cluster_rmsd, best_ENERGY_cluster_rank) ) )

		if(OUTPUT_CHEM_SHIFT): 
			TABLE.write( '%-12s' %(float(get_col_score_value('R_shift_RMSD', best_ENERGY_data_list, best_ENERGY_col_name_list))) )
			TABLE.write( '%-12s' %(float(get_col_score_value('rna_chem_shift', best_ENERGY_data_list, best_ENERGY_col_name_list))) )


		if(DO_output_Rosetta_BP_stats): output_Rosetta_BP_stats( TABLE, best_ENERGY_data_list, best_ENERGY_col_name_list)

		if(DO_output_MC_stats or DO_output_FR3D_stats): 
			output_third_party_stats( TABLE, spec_job, best_ENERGY_data_list, best_ENERGY_col_name_list, cluster_silent_file , "BEST_ENERGY_CLUSTER",  recreate_file)


################################################################################################################################
def extract_best_model(silent_file, col_name): #best here always correspond to the 'lowest'

	print 
	print "extract_best_model FOR: silent_file=(%s) | col_name=(%s)" %(silent_file, col_name)
 

	assert_is_valid_non_empty_silent_file(silent_file)

	filter_outfile="BEST_%s.out" %(col_name)

	###First find the silent_struct tag corresponding to the best model.
	filter_command ="SWA_filter_outfile_wrapper.py -infile %s -filter_outfile %s "  %(silent_file, filter_outfile)
	filter_command+="-max_n_struct 1 -scorecol_name %s -remove_SCORE_file True " %(col_name)
	filter_command+=" > LOG_filter_%s_%s.txt" %(basename(silent_file), col_name)

	submit_subprocess(filter_command)

	assert_is_valid_non_empty_silent_file(filter_outfile)

	###Find the tag of this best model! [also confirm that there is only 1 model.
	data=safe_open(filter_outfile, mode='r' ,Is_master=False) 

	SEQUENCE_LINE = data.readline() # The SEQUENCE: gggcgcagccu line
	COLUMN_NAME_LINE   = data.readline() # The column name line
	col_name_list=COLUMN_NAME_LINE.split()

	struct_count=0
	struct_tag=""

	for line in data:
	
		data_list=line.split()

		if(data_list[0]!='SCORE:'): continue

		struct_count+=1

		struct_tag=data_list[get_col_index_value('description', data_list, col_name_list)]

	if(struct_count!=1): error_exit_with_message("struct_count=(%s)!=1" %(struct_count))
	#############################################################################################
	output_foldername="EXTRACT_%s" %(filter_outfile)

	extract_command ="SWA_extract_pdb.py -silent_file %s -output_foldername %s -tag %s " %(filter_outfile, output_foldername, struct_tag)
	extract_command+=" > LOG_SWA_extract_pdb_%s_%s.txt" %(basename(silent_file), col_name)

	submit_subprocess(extract_command)

	best_pdb="%s.pdb" %(struct_tag)

	if(exists(output_foldername +'/'+ best_pdb)==False):
		error_exit_with_message("%s doesn't exist!" %(output_foldername +'/'+ best_pdb))

	if(exists(best_pdb)==False): #OK possible that this tag was already extracted [since for example best_score_model can be the same struct as best_shift_RMSD_model]
		submit_subprocess("cp %s %s" %(output_foldername +'/'+ best_pdb, best_pdb)) 

	return best_pdb

################################################################################################################################

def get_shift_col_names(silent_file_type):

	shift_RMSD_col_name=""
	shift_SCORE_col_name=""
	num_shift_data_col_name=""

	if(silent_file_type=="MINIMIZED"):
		shift_RMSD_col_name="R_shift_RMSD"
		shift_SCORE_col_name="rna_chem_shift"
		num_shift_data_col_name="R_num_shift_data"
	elif(silent_file_type=="PRE_MINIMIZE"):
		shift_RMSD_col_name="shift_RMSD"
		shift_SCORE_col_name="shift_score"
		num_shift_data_col_name="num_shift_data"
	else:
		error_exit_with_message("Invalid silent_file_type (%s)" %(silent_file_type ))

	return (shift_RMSD_col_name, shift_SCORE_col_name, num_shift_data_col_name)

################################################################################################################################
def output_chem_shift_data_mode(TABLE, col_type, silent_file_type, silent_file_list, common_args_file):

	#col_type_list=['best_ENERGY', 'best_shift_RMSD']	


	(shift_RMSD_col_name, shift_SCORE_col_name, num_shift_data_col_name) = get_shift_col_names(silent_file_type)

	best_model_data_list=[]
	best_col_type_value=99999.99

	rmsd_to_best_shift_col_name='%s_TO_BEST_SHIFT' %(RMSD_COL_NAME_PREFIX)
	rmsd_to_best_energy_col_name='%s_TO_BEST_ENERGY' %(RMSD_COL_NAME_PREFIX)

	col_type_col_name=""
	gap_rmsd_col_name=""

	if(col_type=='best_ENERGY'):
		col_type_col_name='score'
		gap_rmsd_col_name=rmsd_to_best_energy_col_name
	elif(col_type=='best_shift_RMSD'):
		col_type_col_name=shift_RMSD_col_name
		gap_rmsd_col_name=rmsd_to_best_shift_col_name
	else:
		error_exit_with_message("invalid score_col_name (%s)" %(score_col_name))		

	######################################################

	first_cluster_best_energy=99999.99
	first_cluster_best_shift_RMSD=99999.99
	first_cluster_best_rosetta_energy=99999.99

	other_cluster_best_energy=99999.99
	other_cluster_best_shift_RMSD=99999.99
	other_cluster_best_rosetta_energy=99999.99
	######################################################

	for n in range(len(silent_file_list)):

		silent_file=silent_file_list[n]

		assert_is_valid_non_empty_silent_file(silent_file)

		###Find the tag of this best model! [also confirm that there is only 1 model.
		data=safe_open(silent_file, mode='r' ,Is_master=False) 

		SEQUENCE_LINE = data.readline() # The SEQUENCE: gggcgcagccu line
		COLUMN_NAME_LINE   = data.readline() # The column name line
		col_name_list=COLUMN_NAME_LINE.split()

		if(n==0):
			REP_SEQUENCE_LINE=SEQUENCE_LINE
			REP_COLUMN_NAME_LINE=COLUMN_NAME_LINE
			REP_col_name_list=col_name_list
		else:
			if(REP_SEQUENCE_LINE!=SEQUENCE_LINE): error_exit_with_message("REP_SEQUENCE_LINE!=SEQUENCE_LINE")
			if(REP_COLUMN_NAME_LINE!=COLUMN_NAME_LINE): error_exit_with_message("REP_COLUMN_NAME_LINE!=COLUMN_NAME_LINE")
			if(Is_equivalent_list(REP_col_name_list, col_name_list)==False): error_exit_with_message("REP_col_name_list and col_name_list are not equivalent list!")

		struct_count=0
		struct_tag=""

		for line in data:

			data_list=line.split()

			if(data_list[0]!='SCORE:'): continue

			col_type_value=float(get_col_score_value(col_type_col_name, data_list, col_name_list))

			gap_rmsd=float(get_col_score_value(gap_rmsd_col_name, data_list, col_name_list))

			shift_rmsd=float(get_col_score_value(shift_RMSD_col_name, data_list, col_name_list))

			energy=float(get_col_score_value('score', data_list, col_name_list))

			shift_score=float(get_col_score_value(shift_SCORE_col_name, data_list, col_name_list))

			rosetta_energy=energy-shift_score

			if((col_type_value<best_col_type_value)):
				best_col_type_value=col_type_value
				best_model_data_list=copy.deepcopy(data_list)

			if(gap_rmsd<CONFIDENCE_GAP_RMSD):
				if(energy<first_cluster_best_energy): first_cluster_best_energy=energy
				if(rosetta_energy<first_cluster_best_rosetta_energy): first_cluster_best_rosetta_energy=rosetta_energy
				if(shift_rmsd<first_cluster_best_shift_RMSD): first_cluster_best_shift_RMSD=shift_rmsd
			else:
				if(energy<other_cluster_best_energy): other_cluster_best_energy=energy
				if(rosetta_energy<other_cluster_best_rosetta_energy): other_cluster_best_rosetta_energy=rosetta_energy
				if(shift_rmsd<other_cluster_best_shift_RMSD): other_cluster_best_shift_RMSD=shift_rmsd
	
	##########################################################################
	common_args_string=safe_readlines( common_args_file )[0]

	global_sample_res_list=get_option_name_args_safe(common_args_string, "global_sample_res_list", [0])

	fasta_file=get_option_name_args_safe(common_args_string, "fasta", "")

	if(len(global_sample_res_list)==0): error_exit_with_message("global_sample_res_list option does not exist in common_args_string (%s)" %(common_args_string))

	if(fasta_file==""): error_exit_with_message("fasta option does not exist in common_rags (%s)" %(common_args_string))

	sequence = open( fasta_file  ).readlines()[1][:-1]	

	total_res=len(sequence)

	chemical_shift_res_list=add_boundary_seq_num_to_list(global_sample_res_list, total_res, boundary_size=1) 

	##########################################################################

	shift_RMSD_gap_cluster=first_cluster_best_shift_RMSD-other_cluster_best_shift_RMSD

	rosetta_energy_gap_cluster="%6.3f" %((first_cluster_best_rosetta_energy-other_cluster_best_rosetta_energy)/len(chemical_shift_res_list))

	energy_gap_cluster="%6.3f" %((first_cluster_best_energy-other_cluster_best_energy)/len(chemical_shift_res_list))


	###########################################################################
	best_model_shift_RMSD =float(get_col_score_value(shift_RMSD_col_name  , best_model_data_list, REP_col_name_list) )

	best_model_shift_score=float(get_col_score_value(shift_SCORE_col_name , best_model_data_list, REP_col_name_list) )

	best_model_energy     =float(get_col_score_value('score'              , best_model_data_list, REP_col_name_list) ) 

	best_model_rosetta_energy=best_model_energy-best_model_shift_score

	###########################################################################

	shift_RMSD_gap_model=best_model_shift_RMSD-other_cluster_best_shift_RMSD

	rosetta_energy_gap_model="%6.3f" %((best_model_rosetta_energy-other_cluster_best_rosetta_energy)/len(chemical_shift_res_list))

	energy_gap_model="%6.3f" %((best_model_energy-other_cluster_best_energy)/len(chemical_shift_res_list))


	###########################################################################

	#subcol_type_list=["| energy", "rms_NAT", "rms_Bshift", "rms_Benergy", "shift_RMSD", "shift_gapC", "Ros_E_gapC", "E_gapC"]

			
	for subcol_type in subcol_type_list:

		if(subcol_type=='| energy'): 			TABLE.write( '| %-10s' %(best_model_energy) )
		if(subcol_type=='rms_NAT'): 			TABLE.write( '%-12s' %(float(get_col_score_value(RMSD_TO_NAT_COL_NAME					, 	best_model_data_list, REP_col_name_list) ) ) ) 
		if(subcol_type=='rms_Bshift')	: 		TABLE.write( '%-12s' %(float(get_col_score_value(rmsd_to_best_shift_col_name		, best_model_data_list, REP_col_name_list) ) ) ) 
		if(subcol_type=='rms_Benergy'):		TABLE.write( '%-12s' %(float(get_col_score_value(rmsd_to_best_energy_col_name	, best_model_data_list, REP_col_name_list) ) ) ) 
		if(subcol_type=='shift_RMSD'):		TABLE.write( '%-12s' %(best_model_shift_RMSD ) ) 

		if(subcol_type=='shift_gapM'):		TABLE.write( '%-12s' %(100*shift_RMSD_gap_model) ) 
		if(subcol_type=='shift_gapC'):		TABLE.write( '%-12s' %(100*shift_RMSD_gap_cluster) ) 
		if(subcol_type=='Ros_E_gapM'):		TABLE.write( '%-12s' %(rosetta_energy_gap_model) ) 
		if(subcol_type=='Ros_E_gapC'):		TABLE.write( '%-12s' %(rosetta_energy_gap_cluster) ) 
		if(subcol_type=='E_gapM'):				TABLE.write( '%-12s'	%(energy_gap_model) ) 
		if(subcol_type=='E_gapC'):				TABLE.write( '%-12s' %(energy_gap_cluster) ) 


################################################################################################################################

def output_best_score_or_best_rmsd(TABLE, spec_job, col_type, recreate_file):

	have_updated_once=False

	global_rmsd=0.0
	global_score=0.0
	global_silent_file=""

	global_data_list=[]
	global_col_name_list=[]

	for spec_silent_file in spec_job['silent_file_list']:

		lines = open(spec_silent_file).readlines()

		col_name_list=lines[1].split()	

		for line_ID in range(2,len(lines)):

			data_list=lines[line_ID].split()

			if(data_list[0]!='SCORE:'): continue

			#print "col_name_list= ", col_name_list
			#print "data_list= ", data_list

			rmsd=float(get_col_score_value(RMSD_TO_NAT_COL_NAME, data_list, col_name_list))
			score=float(get_col_score_value('score', data_list, col_name_list))

			tag=data_list[get_col_index_value('description', data_list, col_name_list)]
			
			update_best_score_or_rmsd=False

			if(have_updated_once==False): 
				update_best_score_or_rmsd=True

			if(col_type=="best_score"):

				if(score<global_score): 
					update_best_score_or_rmsd=True

			else: #col_type==lowest_rmsd

				if(rmsd<global_rmsd):	 
					update_best_score_or_rmsd=True

			'''
			if(update_best_score_or_rmsd):

				if(IGNORE_FARFAR_NO_AUTO_BULGE_PARENT_TAG):
					PARENT_TAG_LINE=lines[line_ID+1]

					if(PARENT_TAG_LINE[0:18]!="REMARK PARENT_TAG "): error_exit_with_message("PARENT_TAG_LINE[0:18]!=\"REMARK PARENT_TAG \"")

					parent_tag=PARENT_TAG_LINE.split()[2]

					if(parent_tag.count("_NO_AUTO_BULGE")>0):
						WITH_AUTO_BULGE_parent_tag=parent_tag.replace("_NO_AUTO_BULGE", "_WITH_AUTO_BULGE")
						grep_popen_command="grep \"REMARK PARENT_TAG %s\" %s" %(WITH_AUTO_BULGE_parent_tag, spec_silent_file) 	
						grep_lines=popen_and_readlines(grep_popen_command, Is_master=False, tag=("Is_FARFAR_NO_AUTO_BULGE_%s_%s.txt" %(line_ID,basename(spec_silent_file))) )

						if(len(grep_lines)>0):
							if(len(grep_lines)!=1): error_exit_with_message("len(grep_lines)!=1")
							if(grep_lines[0].split()[2]!=WITH_AUTO_BULGE_parent_tag): 
								error_exit_with_message("grep_lines[0].split()[2]=(%s)!=(%s)=WITH_AUTO_BULGE_parent_tag" %(grep_lines[0].split()[2],WITH_AUTO_BULGE_parent_tag))

							print "Ignoring tag (%s) with NO_AUTO_BULGE_parent_tag (%s) since WITH_AUTO_BULGE parent_tag (%s) exist!" %(tag, parent_tag, WITH_AUTO_BULGE_parent_tag)
							update_best_score_or_rmsd=False
			'''

			if(update_best_score_or_rmsd):	

				have_updated_once=True
				global_rmsd=rmsd
				global_score=score
				global_silent_file=spec_silent_file
				global_data_list=copy.deepcopy(data_list)
				global_col_name_list=copy.deepcopy(col_name_list)

	################################		

	TABLE.write( '%-12s' %(global_score) )
	TABLE.write( '%-12s' %(global_rmsd) )


	if(OUTPUT_CHEM_SHIFT): 
		TABLE.write( '%-12s' %(float(get_col_score_value('R_shift_RMSD', global_data_list, global_col_name_list))) )
		TABLE.write( '%-12s' %(float(get_col_score_value('rna_chem_shift', global_data_list, global_col_name_list))) )


	if(DO_output_Rosetta_BP_stats): output_Rosetta_BP_stats( TABLE, global_data_list, global_col_name_list)

	if(DO_output_MC_stats or DO_output_FR3D_stats): 
		output_third_party_stats( TABLE, spec_job, global_data_list, global_col_name_list, global_silent_file, ( "PRE_CLUSTER_%s" %(col_type) ).upper() , recreate_file)


################################################################################################################################


def get_recalculated_silent_file_list(data_src, recal_rmsd_foldername, allow_missing_data): #Update on March 10, 2012

	silent_file_list=["%s/FINAL_calc_RMSD_silent_file.out" %(recal_rmsd_foldername) ]

	silent_file_list.sort()

	if(len(silent_file_list)==0): error_exit_with_message("len(silent_file_list)==0")

	if(allow_missing_data):
		ACT_silent_file_list=[]
		for silent_file in silent_file_list:
			if(exists(silent_file)==False): continue
			ACT_silent_file_list.append(silent_file)

	else:
		ACT_silent_file_list=silent_file_list

	print "------------------------------------------------------------------------------------------"
	print "recalculated_silent_file_list= ", ACT_silent_file_list

	print "------------------------------------------------------------------------------------------"

	for silent_file in ACT_silent_file_list:
		assert_is_valid_non_empty_silent_file(silent_file)

	return ACT_silent_file_list


##############################################################################################################

def get_PRE_silent_file_list(data_src): #Updated on March 10, 2012


	if(data_src['silent_file_type']=="MINIMIZED"):
		#FULL_LENGTH_MINIMIZER/minimize_WITH_SHIFT_STATS_Feb_06_SWA.out

		silent_file_src_folder=get_silent_file_src_foldername(data_src)

		if(exists(silent_file_src_folder)==False):
			error_exit_with_message("silent_file_src_folder (%s) doesn't exist!" %(silent_file_src_folder) )

		PRE_silent_file_list=['%s/FULL_LENGTH_MINIMIZER/minimize_%s' %(silent_file_src_folder, data_src['PRE_minimize_silent_file'])] 

	elif(data_src['silent_file_type']=="PRE_MINIMIZE"):

		PRE_minimize_silent_file= "%s/%s/%s" %(data_src['main_folder'], data_src['folder_name'], data_src['PRE_minimize_silent_file'])
		PRE_silent_file_list=[PRE_minimize_silent_file]

	else:
		error_exit_with_message("Invalid job['silent_file_type']=(%s)" %(data_src['silent_file_type']))

	print 
	print "get_PRE_silent_file_list | silent_file_type=%s | PRE_silent_file_list=%s " %(data_src['silent_file_type'],  list_to_string(PRE_silent_file_list) )

	if(len(PRE_silent_file_list)==0): error_exit_with_message("len(PRE_silent_file_list)==0")

	for n in range(len(PRE_silent_file_list) ): 
		assert_is_valid_non_empty_silent_file(PRE_silent_file_list[n])

	return PRE_silent_file_list

##############################################################################################################

def get_silent_file_src_foldername(data_src): #Written on March 10, 2012

	silent_file_src_foldername=""

	if(data_src["silent_file_folder"][-6:]=="$SELF$"):
		silent_file_src_foldername="%s/" %( data_src["main_folder"][:-6])
	else:
		if(data_src["silent_file_type"]=="MINIMIZED"):
			silent_file_src_foldername="%s/main_folder/%s/%s/" %(  data_src["silent_file_folder"], data_src['algorithm'],  data_src["folder_name"] )
		else:
			silent_file_src_foldername="%s/%s/" %(  data_src["silent_file_folder"],  data_src["folder_name"] )

	if(exists(silent_file_src_foldername)==False): error_exit_with_message("silent_file_src_foldername (%s) doesn't exist " %(silent_file_src_foldername) )

	return silent_file_src_foldername

##############################################################################################################

def find_matching_data_src_list(data_src_list, motif, algorithm): #Written on March 10, 2012

	match_data_src_list=[]

	for curr_data_src in data_src_list:

		if(curr_data_src['folder_name']!=motif): continue

		if((curr_data_src['algorithm']==algorithm) or (algorithm=="SWA_FARFAR_combine")): 
			match_data_src_list.append(copy.deepcopy(curr_data_src))

	return match_data_src_list

##############################################################################################################

def replace_silent_scoreline_func(keep_cols_list, tag_prestring, infile, outfile):

	print
	print "replace_silent_scoreline FOR: infile=(%s) | keep_cols_list=(%s) " %(basename(infile), list_to_string(keep_cols_list, ",")[1:])


	assert_is_valid_non_empty_silent_file(infile)
	if(exists(outfile)): error_exit_with_message("outfile (%s) already  exist!" %(outfile))



	command ="replace_silent_scoreline.py "
	command+=" -keep_column_names %s "%(list_to_string(keep_cols_list))
	command+=" -infile %s " %(infile)
	command+=" -outfile %s " %(outfile)
	if(tag_prestring!=""): command+=" -tag_prestring %s " %(tag_prestring)
	command+=" > LOG_replace_silent_scoreline_%s.txt" %(basename(infile))



	submit_subprocess(command)

	assert_is_valid_non_empty_silent_file(outfile)

##############################################################################################################
def get_working_foldername(silent_file_type, script_algorithm, algorithm, motif):

	working_foldername="%s/%s/%s/%s" %(silent_file_type, script_algorithm.upper(), algorithm, motif)

	working_foldername=os.path.abspath(working_foldername)

	return working_foldername



































