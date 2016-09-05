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

from SWA_dagman_python.parser.SWA_parse_internal_arguments import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
from SWA_dagman_python.parser.SWA_parse_benchmark import parse_benchmark_table_data_file
######################################################################


###Global variables######
print "----------------------------------GLOBAL VARIABLES:----------------------------------"

NO_GRAPHIC=True
perform_VDW_rep_screen=True

TOP_ENERGY_CLUSTERS_FOLDER="TOP_ENERGY_CLUSTERS/"

OUTPUT_TOP_ENERGY_CLUSTERS_MODE=True
GET_OPTIMIZED_NATIVE_MODE=False
OPTIMIZED_NATIVE_RMSD_CUTOFF=1.5 #9999.99 #classify conformation as optimize_native only if RMSD<1.5

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
CHEMICAL_SHIFT_MODE= True

CHEMICAL_SHIFT_SCORE_WEIGHT=4.0

###########################################################

if(OUTPUT_TOP_ENERGY_CLUSTERS_MODE):
	col_type_list=[]	
	for n in range(NUM_ENERGY_CLUSTERS): col_type_list.append("cluster_%d" %(n+1) )
else:
	#col_type_list=["best_score", "lowest_rmsd", "best_RMSD_cluster", "best_ENERGY_cluster" ]
	#Also note implementation below REQUIRES THAT "best_RMSD_cluster" COME RIGHT BEFORE "best_ENERGY_cluster" (BUT OK it either or both is missing)! 
	#col_type_list=["lowest_rmsd", "best_ENERGY_cluster" ]
	col_type_list=["best_score", "lowest_rmsd"]



if(GET_OPTIMIZED_NATIVE_MODE): col_type_list=["best_score"] #don't really need the other two column.



DO_output_Rosetta_BP_stats=False
DO_output_MC_stats=False
DO_output_FR3D_stats=False

subcol_type_list=["energy", "rmsd"]

###########The order the subcol_type are add does matter!#################
if(CHEMICAL_SHIFT_MODE): 
	subcol_type_list.append("shift_RMSD")
	subcol_type_list.append("shift_score")

if(DO_output_Rosetta_BP_stats): subcol_type_list.append("Ros_Pair")
if(DO_output_MC_stats):         subcol_type_list.extend(["MC_Pair","MC_Stack"]) 
if(DO_output_FR3D_stats): 			subcol_type_list.extend(["FR3D_Pair","FR3D_Stack"]) 

INDIV_COL_WIDTH=12
SUB_COLUMNS_WIDTH=len(subcol_type_list)*INDIV_COL_WIDTH
COLUMNS_WIDTH=len(col_type_list)*SUB_COLUMNS_WIDTH

print "--------------------------------------------------------------------"
print "CLUSTER_RMSD= %s " %(CLUSTER_RMSD) 
print "SUITE_CLUSTER_RMSD= %s " %(SUITE_CLUSTER_RMSD)
print "NATIVE_RMSD_CUTOFF= %s " %(NATIVE_RMSD_CUTOFF) 
print "--------------------------------------------------------------------"
print "OUTPUT_TOP_ENERGY_CLUSTERS_MODE=%s" %(GET_OPTIMIZED_NATIVE_MODE)
print "NUM_ENERGY_CLUSTERS= %d " %(NUM_ENERGY_CLUSTERS) 
print "--------------------------------------------------------------------"
print "GET_OPTIMIZED_NATIVE_MODE=%s" %(GET_OPTIMIZED_NATIVE_MODE)
print "OPTIMIZED_NATIVE_RMSD_CUTOFF= %d " %(OPTIMIZED_NATIVE_RMSD_CUTOFF) 
print "--------------------------------------------------------------------"
print "DOUBLE_COUNT_BP_AND_BS=%s" %(DOUBLE_COUNT_BP_AND_BS)
print "INCLUDE_BIFURCARCATION_BP=%s" %(INCLUDE_BIFURCARCATION_BP)
print "--------------------------------------------------------------------"
print "IGNORE_FARFAR_NO_AUTO_BULGE_PARENT_TAG=%s" %(IGNORE_FARFAR_NO_AUTO_BULGE_PARENT_TAG)
print "IGNORE_UNMATCHED_VIRTUAL_RES=%s" %(IGNORE_UNMATCHED_VIRTUAL_RES)
print "--------------------------------------------------------------------"
print "CHEMICAL_SHIFT_MODE=%s " %(CHEMICAL_SHIFT_MODE)
print "CHEMICAL_SHIFT_SCORE_WEIGHT=%s " %(CHEMICAL_SHIFT_SCORE_WEIGHT)
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
def recalculate_silent_file(input_silent_file_list, job, common_args_file, recal_rmsd_silent_file):

	recal_rmsd_command="SWA_cluster.py  -num_pose_kept 9999999 -skip_clustering True -recreate_silent_struct True -extract_pdb False "

	recal_rmsd_command+="-write_score_only False "

	recal_rmsd_command+="-silent_file %s "%(list_to_string(input_silent_file_list))
	recal_rmsd_command+="-native_pdb %s "%(job["native_pdb"])
	recal_rmsd_command+="-common_args %s "%(common_args_file)
	recal_rmsd_command+="-output_filename %s " %(recal_rmsd_silent_file)
	recal_rmsd_command+="-perform_VDW_rep_screen %s " %(perform_VDW_rep_screen)
	recal_rmsd_command+="-no_graphic %s " %(NO_GRAPHIC)
	
	print
	print recal_rmsd_command
	print
	submit_subprocess(recal_rmsd_command)

##############################################################################################################################
def	check_PRE_silent_file(PRE_silent_file):

	if(exists(PRE_silent_file)==False): error_exit_with_message("PRE_silent_file (%s) doesn't exist!: " %(PRE_silent_file) ) 

	if(is_non_empty_silent_file(PRE_silent_file)==False): 
		print "PRE_silent_file (%s) is empty!, moving to next PRE_silent_file!" %(PRE_silent_file) 
		return True

	return False

##############################################################################################################################
def check_for_existing_recalculate_silent_file(recreate_file, recal_rmsd_silent_file):

	if(recreate_file and exists(recal_rmsd_silent_file)==True):
		print "recreate_file==True, recal_rmsd_silent_file (%s) already exist... removing..." %(os.path.abspath(recal_rmsd_silent_file))
		submit_subprocess('rm %s' %(recal_rmsd_silent_file))		

	if(exists(recal_rmsd_silent_file)):
		print "\nrecal_rmsd_silent_file (%s) already exist!\n" %(recal_rmsd_silent_file)
		return True

	return False

##############################################################################################################################
def check_file_list_exist(file_list):
	
	for filename in file_list:
		if(exists(filename)==False): error_exit_with_message("filename (%s) doesn't exist!" %(os.path.abspath(filename)))


##################################################################################################################################

def create_cluster_silent_file(spec_job, good_energy_silent_file_list, recreate_file):

	common_args_file=spec_job["common_args"]

	cluster_silent_file="%s/top_energy_clusters.out" %(TOP_ENERGY_CLUSTERS_FOLDER)

	if(recreate_file and exists(cluster_silent_file)==True):
		print "recreate_file==True, cluster_silent_file (%s) already exist... removing..." %(os.path.abspath(cluster_silent_file))
		submit_subprocess('rm %s' %(cluster_silent_file))		

	if(exists(cluster_silent_file)==False):
		submit_subprocess( "mkdir -p %s " %(TOP_ENERGY_CLUSTERS_FOLDER ) )

		#-write_score_only True only work in recalculate RMSD mode...

		cluster_top_energy_command="SWA_cluster.py  -num_pose_kept 100 -distinguish_pucker false -extract_pdb False " 
		cluster_top_energy_command+="-cluster_rmsd %s " %(CLUSTER_RMSD)
		cluster_top_energy_command+="-suite_cluster_rmsd %s " %(SUITE_CLUSTER_RMSD)
		cluster_top_energy_command+="-silent_file %s "%(list_to_string(good_energy_silent_file_list))
		cluster_top_energy_command+="-native_pdb %s "%(spec_job["native_pdb"])
		cluster_top_energy_command+="-common_args %s "%(common_args_file)
		cluster_top_energy_command+="-output_filename %s " %(cluster_silent_file)	
		cluster_top_energy_command+="-no_graphic %s " %(NO_GRAPHIC)
		cluster_top_energy_command+="-full_length_loop_rmsd_clustering True "

		if(IGNORE_FARFAR_NO_AUTO_BULGE_PARENT_TAG):
			cluster_top_energy_command+="-ignore_FARFAR_no_auto_bulge_parent_tag True "

		if(IGNORE_UNMATCHED_VIRTUAL_RES):
			cluster_top_energy_command+="-ignore_unmatched_virtual_res True "

	########################################################################################################
		submit_subprocess(cluster_top_energy_command)
	else:
		print
		print "cluster_silent_file (%s) already exist!" %(cluster_silent_file)
		print

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

		rmsd=float(get_col_score_value(spec_job['rmsd_col_name'], data_list, col_name_list))
		score=float(get_col_score_value('score', data_list, col_name_list))

		TABLE.write( '%-12s' %(score) )
		TABLE.write( '%-12s' %(rmsd)  )

		if(CHEMICAL_SHIFT_MODE): 
			TABLE.write( '%-12s' %(float(get_col_score_value('shift_RMSD', data_list, col_name_list))) )
			TABLE.write( '%-12s' %(float(get_col_score_value('shift_score', data_list, col_name_list))) )


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

		rmsd=float(get_col_score_value(spec_job['rmsd_col_name'], data_list, col_name_list))
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

		if(CHEMICAL_SHIFT_MODE): 
			TABLE.write( '%-12s' %(float(get_col_score_value('shift_RMSD', best_RMSD_data_list, best_RMSD_col_name_list))) )
			TABLE.write( '%-12s' %(float(get_col_score_value('shift_score', best_RMSD_data_list, best_RMSD_col_name_list))) )


		if(DO_output_Rosetta_BP_stats): output_Rosetta_BP_stats( TABLE, best_RMSD_data_list, best_RMSD_col_name_list)

		if(DO_output_MC_stats or DO_output_FR3D_stats): 
			output_third_party_stats( TABLE, spec_job, best_RMSD_data_list, best_RMSD_col_name_list, cluster_silent_file, "BEST_RMSD_CLUSTER",  recreate_file)

	if(col_type_list.count("best_ENERGY_cluster")==1):
		TABLE.write( '%-12s' %(best_ENERGY_cluster_score) )
		if(not_enough_cluster):
			TABLE.write( '%-12s' %('%s(!%d)' %(best_ENERGY_cluster_rmsd, best_ENERGY_cluster_rank) ) )
		else:
			TABLE.write( '%-12s' %('%s(%d)' %(best_ENERGY_cluster_rmsd, best_ENERGY_cluster_rank) ) )

		if(CHEMICAL_SHIFT_MODE): 
			TABLE.write( '%-12s' %(float(get_col_score_value('shift_RMSD', best_ENERGY_data_list, best_ENERGY_col_name_list))) )
			TABLE.write( '%-12s' %(float(get_col_score_value('shift_score', best_ENERGY_data_list, best_ENERGY_col_name_list))) )


		if(DO_output_Rosetta_BP_stats): output_Rosetta_BP_stats( TABLE, best_ENERGY_data_list, best_ENERGY_col_name_list)

		if(DO_output_MC_stats or DO_output_FR3D_stats): 
			output_third_party_stats( TABLE, spec_job, best_ENERGY_data_list, best_ENERGY_col_name_list, cluster_silent_file , "BEST_ENERGY_CLUSTER",  recreate_file)


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

			rmsd=float(get_col_score_value(spec_job['rmsd_col_name'], data_list, col_name_list))
			score=float(get_col_score_value('score', data_list, col_name_list))

			tag=data_list[get_col_index_value('description', data_list, col_name_list)]

			if(GET_OPTIMIZED_NATIVE_MODE):
				if(rmsd>=OPTIMIZED_NATIVE_RMSD_CUTOFF): continue

			
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


	if(CHEMICAL_SHIFT_MODE): 
		TABLE.write( '%-12s' %(float(get_col_score_value('shift_RMSD', global_data_list, global_col_name_list))) )
		TABLE.write( '%-12s' %(float(get_col_score_value('shift_score', global_data_list, global_col_name_list))) )


	if(DO_output_Rosetta_BP_stats): output_Rosetta_BP_stats( TABLE, global_data_list, global_col_name_list)

	if(DO_output_MC_stats or DO_output_FR3D_stats): 
		output_third_party_stats( TABLE, spec_job, global_data_list, global_col_name_list, global_silent_file, ( "PRE_CLUSTER_%s" %(col_type) ).upper() , recreate_file)


################################################################################################################################


def get_recalculated_silent_file_list(spec_job, GET_OPTIMIZED_NATIVE_MODE):

	good_energy_silent_file_list=[]

	WITH_chem_shift_str=""
	if(CHEMICAL_SHIFT_MODE): WITH_chem_shift_str="WITH_SHIFT_STATS_"

	if(spec_job['silent_file']=="QUICK"):

		if( spec_job['algorithm'].count("FARFAR")!=0 ): 
			error_exit_with_message('"spec_job[\'silent_file\']=="QUICK" but job[\'algorithm\'].count("FARFAR")!=0" ' )


		spec_job['silent_file_list'] = glob( "FILTERED_RECAL_RMSD/recal_filtered_score*.out")

		good_energy_silent_file_list=copy.deepcopy(spec_job['silent_file_list'])

		if(GET_OPTIMIZED_NATIVE_MODE==False): spec_job['silent_file_list'].append("FILTERED_RECAL_RMSD/recal_filtered_RMSD*.out")

	elif(spec_job['silent_file']=="QUICK_INDIVIDUAL_REGIONS"):

		if( spec_job['algorithm'].count("FARFAR")!=0 ): 
			error_exit_with_message('"spec_job[\'silent_file\']=="QUICK_INDIVIDUAL_REGIONS" but job[\'algorithm\'].count("FARFAR")!=0" ' )

		spec_job['silent_file_list']=["INDIVIDUAL_REGIONS/FILTERED_RECAL_RMSD/recal_filtered_score_individual_regions.out"]

		good_energy_silent_file_list=copy.deepcopy(spec_job['silent_file_list'])

		if(GET_OPTIMIZED_NATIVE_MODE==False): spec_job['silent_file_list'].append("INDIVIDUAL_REGIONS/FILTERED_RECAL_RMSD/recal_filtered_RMSD_individual_regions.out")

	elif(spec_job['silent_file']=="STANDARD"):

		if( spec_job['algorithm'].count("FARFAR")!=0 ): 

			recal_FARFAR_score='RECAL_RMSD/recal_full_%sFINAL_filtered_energy.out' %(WITH_chem_shift_str)
			recal_FARFAR_RMSD='RECAL_RMSD/recal_full_%sFINAL_filtered_RMSD.out' %(WITH_chem_shift_str)

			'''Comment out on Nov 06, 2011!
			##########Consistency Check#############################################
			glob_test=glob( "RECAL_RMSD/recal_full_*.out" )

			if(len(glob_test)!=2): 
				print "glob_test= ", glob_test 
				error_exit_with_message("len(glob_test)!=2")

			if(glob_test[0]!=recal_FARFAR_score): 
				print "glob_test= ", glob_test 
				error_exit_with_message("glob_test[0]!=recal_FARFAR_score]")

			if(glob_test[1]!=recal_FARFAR_RMSD):  
				print "glob_test= ", glob_test
				error_exit_with_message("glob_test[1]!=recal_FARFAR_RMSD)")

			########################################################################
			'''

			spec_job['silent_file_list'] = [recal_FARFAR_score]
			good_energy_silent_file_list = [recal_FARFAR_score]

			if(GET_OPTIMIZED_NATIVE_MODE==False): spec_job['silent_file_list'].append(recal_FARFAR_RMSD)

		else:

			#spec_job['silent_file_list'] = glob( "RECAL_RMSD/recal_full_*.out" ) Comment out on Nov 06, 2011
			spec_job['silent_file_list']=["RECAL_RMSD/recal_full_%sregion_FINAL.out" %(WITH_chem_shift_str)]
			good_energy_silent_file_list=copy.deepcopy(spec_job['silent_file_list'])

	elif(spec_job['silent_file']=="INDIVIDUAL_REGIONS"):

		if( spec_job['algorithm'].count("FARFAR")!=0 ): 
			error_exit_with_message('"spec_job[\'silent_file\']=="INDIVIDUAL_REGIONS" but job[\'algorithm\'].count("FARFAR")!=0" ' )

		spec_job['silent_file_list']=[]
		full_length_REGION_FOLDER_list=get_full_length_REGION_FOLDER_list("CONDOR/REGION_FINAL_cluster.condor", verbose=True)
		for REGION_FOLDER in full_length_REGION_FOLDER_list:

			globstring = '%s/RECAL_RMSD/recal_full_start_from_region_*_sample_filtered.out' %(REGION_FOLDER)
			spec_job['silent_file_list'].extend( glob( globstring ) )

		good_energy_silent_file_list=copy.deepcopy(spec_job['silent_file_list'])

	elif(spec_job['silent_file'].count("SWA")!=0): #New Nov 06, 2011

		if( spec_job['algorithm'].count("FARFAR")>0 ): 	error_exit_with_message("spec_job['algorithm'].count(\"FARFAR\")>0!") 

		if(spec_job['silent_file']=="SWA_RB_STANDARD"):

			spec_job['silent_file_list']=["FINAL_REBUILD_BULGE/STANDARD/RECAL_RMSD/recal_full_%srebuild_bulge_region_FINAL.out" %(WITH_chem_shift_str)] 

		else:
			error_exit_with_message("Invalid SWA spec_job['silent_file']=%s" %(spec_job['silent_file']))

		good_energy_silent_file_list=copy.deepcopy(spec_job['silent_file_list'])


	elif(spec_job['silent_file'].count("FARFAR")!=0): #New Nov 06, 2011

		if( spec_job['algorithm'].count("FARFAR")==0 ): 	error_exit_with_message("spec_job['algorithm'].count(\"FARFAR\")==0!") 

		recal_FARFAR_score="EMPTY"
		recal_FARFAR_RMSD="EMPTY"

		if(spec_job['silent_file']=="FARFAR_CLUSTERED"):
			recal_FARFAR_score='RECAL_RMSD/recal_full_%sclustered_FINAL_filtered_energy.out' %(WITH_chem_shift_str)
			recal_FARFAR_RMSD='RECAL_RMSD/recal_full_%sclustered_FINAL_filtered_RMSD.out' %(WITH_chem_shift_str)

		elif(spec_job['silent_file']=="FARFAR_RB_STANDARD"):
			FINAL_filter_RMSD_file='FINAL_RMSD_REBUILD_BULGE/STANDARD/RECAL_RMSD/recal_full_%srebuild_bulge_clustered_FINAL_filtered_RMSD.out' %(WITH_chem_shift_str)
			FINAL_filter_energy_file='FINAL_ENERGY_REBUILD_BULGE/STANDARD/RECAL_RMSD/recal_full_%srebuild_bulge_clustered_FINAL_filtered_energy.out' %(WITH_chem_shift_str)

		elif(spec_job['silent_file']=="FARFAR_RB_WITH_SHIFT"):
			FINAL_filter_RMSD_file='FINAL_RMSD_REBUILD_BULGE/WITH_CHEM_SHIFT/RECAL_RMSD/recal_full_%srebuild_bulge_clustered_FINAL_filtered_RMSD.out' %(WITH_chem_shift_str)
			FINAL_filter_energy_file='FINAL_ENERGY_REBUILD_BULGE/WITH_CHEM_SHIFT/RECAL_RMSD/recal_full_%srebuild_bulge_clustered_FINAL_filtered_energy.out' %(WITH_chem_shift_str)


		else:
			error_exit_with_message("Invalid FARFAR spec_job['silent_file']=%s" %(spec_job['silent_file']))

		spec_job['silent_file_list'] = [recal_FARFAR_score]
		good_energy_silent_file_list = [recal_FARFAR_score]

		if(GET_OPTIMIZED_NATIVE_MODE==False): spec_job['silent_file_list'].append(recal_FARFAR_RMSD)


	else:
		error_exit_with_message("Invalid spec_job['silent_file'] (%s) " %(spec_job['silent_file']) )

	spec_job['silent_file_list'].sort()
	good_energy_silent_file_list.sort()

	check_file_list_exist(copy.deepcopy(spec_job['silent_file_list']))
	check_file_list_exist(copy.deepcopy(good_energy_silent_file_list))


	print "------------------------------------------------------------------------------------------"
	print "spec_job['silent_file_list']= ", spec_job['silent_file_list']
	print "good_energy_silent_file_list= ", good_energy_silent_file_list
	print "------------------------------------------------------------------------------------------"

	#if(len(spec_job['silent_file_list'])==0): error_exit_with_message("len(spec_job['silent_file_list'])==0")
	#if(len(good_energy_silent_file_list)==0): error_exit_with_message("len(good_energy_silent_file_list)==0")


	return (spec_job['silent_file_list'], good_energy_silent_file_list)


##############################################################################################################

def get_PRE_silent_file_list(job, quick):


	if(quick and (job['algorithm'].count("FARFAR")!=0) ):
		error_exit_with_message('"quick==True but job[\'algorithm\'].count("FARFAR")!=0, job[\'algorithm\']==%s' %(job['algorithm'])  )

	if(job['PRE_silent_file']=="INDIVIDUAL_REGIONS"):
		if( job['algorithm'].count("FARFAR")!=0 ): 
			error_exit_with_message('"job[\'PRE_silent_file\']=="INDIVIDUAL_REGIONS" but job[\'algorithm\'].count("FARFAR")!=0, job[\'algorithm\']==%s' %(job['algorithm']) )

		job['PRE_silent_file_list']=[]

		final_cluster_dag_file="CONDOR/REGION_FINAL_cluster.condor"

		if(exists(final_cluster_dag_file)==False): 
			final_cluster_dag_file="CONDOR/CLUSTERER/REGION_FINAL_cluster.condor"

		if(exists(final_cluster_dag_file)==False): 
			error_exit_with_message("final_cluster_dag_file (%s) doesn't exist!" %(final_cluster_dag_file))

		full_length_REGION_FOLDER_list=get_full_length_REGION_FOLDER_list(final_cluster_dag_file, verbose=False)
		for REGION_FOLDER in full_length_REGION_FOLDER_list:

			globstring = '%s/start_from_region_*_sample_filtered.out' %(REGION_FOLDER) #Used to be start_from_region_*_*_sample_filtered.out until March 28, 2011
			glob_file_list = glob( globstring )
			glob_file_list.sort()
			job['PRE_silent_file_list'].extend(glob_file_list)

	elif(job['PRE_silent_file'].count("FARFAR")>0):

		if( job['algorithm'].count("FARFAR")==0 ): 
			error_exit_with_message('"job[\'PRE_silent_file\']=="FARFAR_DAGMAN" but job[\'algorithm\']=%s ' %(job['algorithm'])  )

		FINAL_filter_RMSD_file=""
		FINAL_filter_energy_file=""

		if(job['PRE_silent_file']=="FARFAR_DAGMAN"):
			FINAL_filter_RMSD_file='FINAL_filtered_RMSD.out'
			FINAL_filter_energy_file='FINAL_filtered_energy.out'

		if(job['PRE_silent_file']=="FARFAR_CLUSTERED"):
			FINAL_filter_RMSD_file='clustered_FINAL_filtered_RMSD.out'
			FINAL_filter_energy_file='clustered_FINAL_filtered_energy.out'

		elif(job['PRE_silent_file']=="FARFAR_RB_STANDARD"):
			FINAL_filter_RMSD_file='FINAL_RMSD_REBUILD_BULGE/STANDARD/rebuild_bulge_clustered_FINAL_filtered_RMSD.out'
			FINAL_filter_energy_file='FINAL_ENERGY_REBUILD_BULGE/STANDARD/rebuild_bulge_clustered_FINAL_filtered_energy.out'

		elif(job['PRE_silent_file']=="FARFAR_RB_WITH_SHIFT"):
			FINAL_filter_RMSD_file='FINAL_RMSD_REBUILD_BULGE/WITH_CHEM_SHIFT/rebuild_bulge_clustered_FINAL_filtered_RMSD.out'
			FINAL_filter_energy_file='FINAL_ENERGY_REBUILD_BULGE/WITH_CHEM_SHIFT/rebuild_bulge_clustered_FINAL_filtered_energy.out'
		else:
			error_exit_with_message("Invalid FARFAR job['PRE_silent_file']=%s" %(job['PRE_silent_file']))

		if(exists(FINAL_filter_RMSD_file)==False): error_exit_with_message("FINAL_filter_RMSD_file (%s) doesn't exist" %(FINAL_filter_RMSD_file)  )

		if(exists(FINAL_filter_energy_file)==False): error_exit_with_message("FINAL_filer_energy_file (%s) doesn't exist" %(FINAL_filter_energy_file) )

		job['PRE_silent_file_list']=[FINAL_filter_RMSD_file, FINAL_filter_energy_file]

	elif(job['PRE_silent_file']=="SWA_RB_STANDARD"):

		if( job['algorithm'].count("FARFAR")>0 ): 	
			print "ERROR job: ", job
			error_exit_with_message("job['algorithm'].count(\"FARFAR\")>0!") 

		job['PRE_silent_file_list']=["FINAL_REBUILD_BULGE/STANDARD/rebuild_bulge_region_FINAL.out"]

	else:
		job['PRE_silent_file_list']=[job['PRE_silent_file']]

	print "------------------------------------------------------------------------------------------"
	print "job['PRE_silent_file']= ", job['PRE_silent_file']
	print "job['PRE_silent_file_list']= ", job['PRE_silent_file_list']
	print "------------------------------------------------------------------------------------------"

	if(len(job['PRE_silent_file_list'])==0): error_exit_with_message("len(job['PRE_silent_file_list'])==0")

	for PRE_silent_file in job["PRE_silent_file_list"]:
		if(exists(PRE_silent_file)==False): error_exit_with_message("PRE_silent_file (%s) doesn't exist!" %(PRE_silent_file))

##############################################################################################################
###Nov 05, 2011
def create_silent_file_WITH_chemical_shift_stats(job, start_silent_file, working_foldername, recal_shifts):

	print "Adding chemical shift stats to start_silent_file (%s)" %(start_silent_file)

	if(exists(working_foldername)==False): error_exit_with_message("working_foldername (%s) doesn't exist!" %(working_foldername))

	base_folder=os.path.abspath(".")

	start_silent_file=os.path.abspath(start_silent_file)

	os.chdir( working_foldername )

	####################################################################################################################
	if(job.has_key("BMRB_file")==False): 
		print "ERROR JOB: ", job
		error_exit_with_message("job.has_key(\"BMRB_file\")==False!")

	BMRB_chemical_shift_file=job["BMRB_file"]
	chemical_shift_args_file=job["src_folder"] + "/" +  get_chemical_shift_args_file()
	common_args_file=job["common_args"]

	if(exists(BMRB_chemical_shift_file)==False): error_exit_with_message("BMRB_chemical_shift_file (%s) doesn't exist!" %(BMRB_chemical_shift_file))

	if(exists(chemical_shift_args_file)==False): error_exit_with_message("chemical_shift_args_file (%s) doesn't exist!")

	if(exists(common_args_file)==False): error_exit_with_message("common_args_file (%s) doesn't exist!" %(common_args_file))

	common_args_list=safe_readlines( common_args_file )[0].split()

	global_sample_res_list=parse_options( common_args_list, "global_sample_res_list", [0])
	if(len(global_sample_res_list)==0): error_exit_with_message("global_sample_res_list option doesn't exist in common_args: %s!" %(list_to_string(common_args)) )

	fasta_file=parse_options( common_args_list, "fasta", "")
	if(fasta_file==""): error_exit_with_message("fasta==\"\"")

	fasta_file=job["src_folder"] + fasta_file

	if(exists(fasta_file)==False): error_exit_with_message("fasta_file (%s) doesn't exist!" %(fasta_file) )

	total_res = len(open( fasta_file ).readlines()[1][:-1])

	calc_chemical_shift_res_list=add_boundary_seq_num_to_list(global_sample_res_list, total_res, boundary_size=1) 

	output_silent_file="WITH_SHIFT_STATS_%s" %(basename(start_silent_file))

	if(exists(output_silent_file)): 
		print "create_silent_file_WITH_chemical_shift_stats::WARNING output_silent_file (%s) already exist!....removing.." %(output_silent_file) 
		submit_subprocess("rm %s " %(output_silent_file))

	add_chemical_shift_stats_command="add_chemical_shift_score_to_silent_file.py "
	add_chemical_shift_stats_command+="-input_BMRB %s " %(BMRB_chemical_shift_file)
	add_chemical_shift_stats_command+="-residue_list %s " %(list_to_string(calc_chemical_shift_res_list))
	add_chemical_shift_stats_command+="-input_silent_file %s " %(start_silent_file)
	add_chemical_shift_stats_command+="-output_silent_file %s " %(output_silent_file)
	add_chemical_shift_stats_command+="-recalculate_shifts %s " %(recal_shifts) 
	add_chemical_shift_stats_command+="-shift_score_weight %s " %(CHEMICAL_SHIFT_SCORE_WEIGHT)
	add_chemical_shift_stats_command+=parse_chemical_shift_args( chemical_shift_args_file )

	#add_chemical_shift_stats_command+=" > LOG_add_shift_stats.out 2> LOG_add_shift_stats.err "

	submit_subprocess(add_chemical_shift_stats_command)
	####################################################################################################################

	os.chdir( base_folder )

	output_silent_file="%s/%s" %(working_foldername, output_silent_file)

	if(exists(output_silent_file)==False): 
		error_exit_with_message("create_silent_file_WITH_chemical_shift_stats::output_silent_file (%s) doesn't exist!" %(output_silent_file))

	return output_silent_file

####################################################################################################################

