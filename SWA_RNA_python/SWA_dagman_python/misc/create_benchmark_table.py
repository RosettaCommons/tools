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

######################################################################

from create_benchmark_table_util import *

from SWA_dagman_python.misc.SWA_cat_outfiles import concatenate_outfiles

from SWA_dagman_python.utility.SWA_util import *

from SWA_dagman_python.parser.SWA_parse_internal_arguments import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options, get_option_name_args_safe
from SWA_dagman_python.parser.SWA_parse_benchmark import *
######################################################################



###SEPT 29, 2011. NOTE TO SELF SWA_filter_outfile_wrapper actually doesn't generate exactly "max_n_struct" structures. This is due to equivalency of the energy score are the cutoff (HAVE TO FIX THIS BEFORE the next time that I run the create_benchmark_table.py code!)

#create_benchmark_table.py -algorithm recal_rmsd -recreate_file False  -star_mode False  -pre_minimize_chem_shift_weight 0.0  > log_recal_rmsd.out 2> log_recal_rmsd.err

#create_benchmark_table.py -algorithm create_table  -recreate_file False -star_mode False  -allow_missing_data True 
#> create_table.out 2> create_table.err



job_submission_file="/Users/sripakpa/minirosetta/Rosetta_rna_file/CS_PAPER/CS_benchmark_rna_list.txt"
data_src_file="/Users/sripakpa/minirosetta/Rosetta_rna_file/CS_PAPER/BENCHMARK_TABLE/DATA_SRC/March_05_2012_CS_benchmark_data_src.txt"
motif_file="/Users/sripakpa/minirosetta/Rosetta_rna_file/CS_PAPER/BENCHMARK_TABLE/DATA_SRC/CS_benchmark_motif_list.txt"

if(exists(job_submission_file)==False): error_exit_with_message("job_submission_file (%s) doesn't exist!!" %(job_submission_file))  
if(exists(data_src_file)==False): error_exit_with_message("data_src_file (%s) doesn't exist!!" %(data_src_file))  
if(exists(motif_file)==False): error_exit_with_message("motif_file (%s) doesn't exist!!" %(motif_file))  

print "job_submission_file=%s" %(job_submission_file)
print "data_src_file=%s" %(data_src_file)
print "motif_file=%s" %(motif_file)


allow_missing_data= parse_options( argv, "allow_missing_data", "False" )
recreate_file= parse_options( argv, "recreate_file", "True" )
script_algorithm= parse_options( argv, "algorithm", "create_table" )
quick= parse_options( argv, "quick", "False" )
star_mode= parse_options( argv, "star_mode", "False" )
pre_minimize_chem_shift_weight= parse_options( argv, "pre_minimize_chem_shift_weight", 4.0 )


if(script_algorithm!="create_table" and script_algorithm!="recal_rmsd"):  error_exit_with_message("Invalid script_algorithm: %s" %(script_algorithm)) 

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )


######################################################


motif_list=parse_motif_file(motif_file)

(data_src_list, algorithm_list)=parse_benchmark_data_source_file(data_src_file, star_mode, parse_silent_file=True)

job_submission_list=parse_benchmark_job_file(job_submission_file, setup_folders_and_files=False, star_mode=False, verbose=False) #ALWAYS set star_mode to false!

print "------------------------------------------------------------------------------------------"
print "------------------------------------------------------------------------------------------"
print "------------------------------------------------------------------------------------------"
print "------------------------------------------------------------------------------------------"

START_algorithm_list=copy.deepcopy(algorithm_list)
algorithm_list=["SWA_FARFAR_combine"]
#algorithm_list.extend(START_algorithm_list)

START_FOLDER=os.path.abspath(".")

if(script_algorithm=="create_table"):

	table_filename="TABLE.txt"

	TABLE = open( table_filename, 'w')
	TABLE.write( '%s\n' %(data_src_file) )

	print_column_names(TABLE, algorithm_list)

	for motif in motif_list:

		(job, found_job)=find_matching_job(job_submission_list, motif)

		if(found_job==False): error_exit_with_message("Cannot fild job for motif=%s " %(motif))

		if(DOUBLE_LINE_SPACING): TABLE.write( '\n' )
		TABLE.write( '%-34s' %(motif) )

		for algorithm in algorithm_list:

			os.chdir(START_FOLDER)

			match_data_src_list=find_matching_data_src_list(data_src_list, motif, algorithm)

			############################################				
			if(len(match_data_src_list)==0): 
				print "motif=%s | algorithm=%s doesn't exist in the data_src_list" %( motif, algorithm)
				for n in range(len(col_type_list)):	TABLE.write( '%-*s' %(COLUMNS_WIDTH,'') )
				continue

			REP_data_src=match_data_src_list[0] #Representative data_src

			silent_file_type=REP_data_src['silent_file_type']
				
			working_foldername=get_working_foldername(silent_file_type, script_algorithm, algorithm, motif)

			recal_rmsd_foldername=get_working_foldername(silent_file_type, "recal_rmsd", algorithm, motif)

			############################################

			print "------------------------------------------------------------------------------------------"
			print "------------------------------------------------------------------------------------------"
			print "extracting data for motif=%s | algorithm=%s | len(match_data_src_list)=%s :" %(motif, algorithm, len(match_data_src_list))
			print 

			if(recreate_file and exists(working_foldername)): #This clean up working_foldername from prior run if recreate_file == True..
				print "working_foldername (%s) already exist ... removing!" %(working_foldername)
				submit_subprocess( "rm -r %s " %(working_foldername) ) #Might as well delete the whole folder since will recreate the file

			submit_subprocess( "mkdir -p %s " %(working_foldername) )
			os.chdir( working_foldername )

			#################################################

			native_pdb=job['PRISTINE_motif_folder_name'] + '/' + job['native_pdb']
			common_args_file= "%s/%s/%s" %(REP_data_src['main_folder'], REP_data_src['folder_name'], REP_data_src['common_args'])

			if(exists(native_pdb)==False): error_exit_with_message("native_pdb (%s) doesn't exist!: " %(native_pdb) ) 
			if(exists(common_args_file)==False): error_exit_with_message("common_args_file (%s) doesn't exist!" %(common_args_file))	

			####common_args_file implicitly need the VDW_rep_screener PDBs and fasta#####
			create_fasta_from_data_src(REP_data_src, native_pdb)

			#if(PERFORM_VDW_REP_SCREEN): create_VDW_rep_pdb_from_data_src(REP_data_src)
			#################################################

			silent_file_list=get_recalculated_silent_file_list(copy.deepcopy(REP_data_src), recal_rmsd_foldername, allow_missing_data)

			if(len(silent_file_list)==0): #Basically this condition is possible only if allow_missing_data==True
				print "len(silent_file_list)==0!"
				for n in range(len(col_type_list)):	TABLE.write( '%-*s' %(COLUMNS_WIDTH,'') )
				continue



			#########Create the cluster_top_energy_silent_file######################################################
			if(CHEM_SHIFT_DATA_MODE):

				for col_type in col_type_list:
					output_chem_shift_data_mode(TABLE, col_type, silent_file_type, silent_file_list, common_args_file)



			else:
				cluster_silent_file=""
				SCORE_cluster_silent_file=""

				if( ( "best_RMSD_cluster" in col_type_list) or ("best_ENERGY_cluster" in col_type_list ) or OUTPUT_TOP_ENERGY_CLUSTERS_MODE):

					cluster_silent_file=create_cluster_silent_file(native_pdb, common_args_file, REP_data_src['silent_file_list'], recreate_file)

					SCORE_cluster_silent_file=append_to_basename("SCORE_", cluster_silent_file)

					submit_subprocess('grep "SCORE: " %s > %s' %(cluster_silent_file, SCORE_cluster_silent_file) )

				########################################################################################################
				already_outputted_best_cluster_stat=False

				for col_type in col_type_list:

					if(col_type=="best_RMSD_cluster" or col_type=="best_ENERGY_cluster" or OUTPUT_TOP_ENERGY_CLUSTERS_MODE): 
					#Need to output this AFTER "best_score" and "lowest_rmsd" to match col_name!! (i.e. follow the col_type_list order!)

						if(already_outputted_best_cluster_stat): continue
						already_outputted_best_cluster_stat=True
					
						if(OUTPUT_TOP_ENERGY_CLUSTERS_MODE):
							output_top_energy_clusters_stats(TABLE, SCORE_cluster_silent_file, cluster_silent_file, REP_data_src, recreate_file)

						else:
							output_best_energy_and_best_RMSD_cluster_stats(TABLE, SCORE_cluster_silent_file, cluster_silent_file, REP_data_src, col_type_list, recreate_file)

					else: #col_type=="best_score" or col_type=="lowest_rmsd"

						output_best_score_or_best_rmsd(TABLE, REP_data_src, col_type, recreate_file)

			TABLE.flush()

		TABLE.write('\n')

	TABLE.close()

else: #recal_rmsd 

	print "------------------------------------------------------------------------------------------"
	print "------------------------------------------------------------------------------------------"
	if(script_algorithm!="recal_rmsd"): error_exit_with_message("Invalid script_algorithm (%s) " %(script_algorithm))

	for algorithm in algorithm_list: #Move algorithm_list to outer loop, so that will finish SWA_FARFAR_combine jobs before moving to the individual SWA_denovo and FARFAR_denovo jobs.

		for motif in motif_list:

			(job, found_job)=find_matching_job(job_submission_list, motif)

			if(found_job==False): error_exit_with_message("Cannot fild job for motif=%s " %(motif))


			os.chdir(START_FOLDER)

			match_data_src_list=find_matching_data_src_list(data_src_list, motif, algorithm)

			############################################				
			if(len(match_data_src_list)==0): 
				print "motif=%s | algorithm=%s doesn't exist in the data_src_list" %( motif, algorithm)
				continue

			REP_data_src=match_data_src_list[0] #REPRESENATIVE data_src.

			silent_file_type=REP_data_src['silent_file_type']

			working_foldername=get_working_foldername(silent_file_type, script_algorithm, algorithm, motif)

			########################################################
			##ASSUME CHEMICAL_SHIFT_DATA BENCHMARK!

			print "------------------------------------------------------------------------------------------"
			print "------------------------------------------------------------------------------------------"
			print "recalculating rmsd for motif=%s | algorithm=%s | len(match_data_src_list)=%s :" %(motif, algorithm, len(match_data_src_list))


			if(recreate_file and exists(working_foldername)): ##IF recreate_file then delete UNDER ALL circumstances
				print "working_foldername (%s) already exist ... removing!" %(working_foldername)
				submit_subprocess( "rm -r %s " %(working_foldername) ) 

			FINAL_CALC_silent_file="FINAL_calc_RMSD_silent_file.out"

			if(exists(working_foldername + '/' + FINAL_CALC_silent_file)): 
				print "FINAL_CALC_silent_file (%s) already exist ... SKIP!" %(FINAL_CALC_silent_file)
				continue
			else: 
				if(exists(working_foldername)):
					##IF FINAL_CALC_silent_file doesn't exist, then INCOMPLETE script run on this motif, NEED TO RERUN FROM SCRATCH!
					print "working_foldername (%s) already exist ... removing!" %(working_foldername)
					submit_subprocess( "rm -r %s " %(working_foldername) ) #Might as well delete the whole folder since will recreate the file


			submit_subprocess( "mkdir -p %s " %(working_foldername) )
			os.chdir( working_foldername )

			TEMP_FOLDER="TEMP/"
			if(exists(TEMP_FOLDER)): error_exit_with_message("TEMP_FOLDER (%s) already exist!")
			submit_subprocess( "mkdir -p %s " %(TEMP_FOLDER) )
			os.chdir( TEMP_FOLDER )

			############################################

			for data_src_ID in range(len(match_data_src_list)):			

				data_src=match_data_src_list[data_src_ID]

				#print "------------------------------------------------------------------------------------------"
				#print "data_src= ", data_src


				data_src['PRE_silent_file_list']=get_PRE_silent_file_list(copy.deepcopy(data_src)) 

				match_data_src_list[data_src_ID]=data_src #Important since 'PRE_silent_file_list' is changed in the for loop!


			#################################################


			native_pdb=job['PRISTINE_motif_folder_name'] + '/' + job['native_pdb']
			common_args_file= "%s/%s/%s" %(REP_data_src['main_folder'], REP_data_src['folder_name'], REP_data_src['common_args'])

			if(exists(native_pdb)==False): error_exit_with_message("native_pdb (%s) doesn't exist!: " %(native_pdb) ) 
			if(exists(common_args_file)==False): error_exit_with_message("common_args_file (%s) doesn't exist!" %(common_args_file))

			####common_args_file implicitly need the VDW_rep_screener PDBs and fasta#####
			create_fasta_from_data_src(REP_data_src, native_pdb)

			if(PERFORM_VDW_REP_SCREEN): create_VDW_rep_pdb_from_data_src(REP_data_src)


			########################################################
			#1. Copy pre_silent_files in working_folder
			#2. Add 1JJ2 or SWA prefix and replace_score_col keep O_loop_rmsd rna_chem_shift R_shift_RMSD R_num_shift_data (and possibly shift_RMSD shift_score num_shift_data)
			#2. CONCATENATE all PRE_silent_files 
			#3. Calculate RMSD to native_pdb. (keep and rename NEW_Full_L_rmsd AND NEW_NAT_rmsd to RMSD_V1_TO_NAT and RMSD_V2_TO_NAT)
			#4. Find model with lowest R_shift_RMSD 
			#5. Calculate RMSD to tho this model (keep and rename NEW_Full_L_rmsd AND NEW_NAT_rmsd to RMSD_V1_TO_BEST_SHIFT and RMSD_V2_TO_BEST_SHIFT)
			#6. Find model with lowest total energy
			#7. Calculate RMSD to tho this model (keep and rename NEW_Full_L_rmsd AND NEW_NAT_rmsd to RMSD_V1_TO_BEST_ENERGY and RMSD_V2_TO_BEST_ENERGY)
			######################################################## 

			#SWA_cat_outfiles_wrapper.py -add_file_num_to_tag False -outfile WITH_SHIFT_STATS_UNIQUE_H5_Feb_10_SWA_Dec_01_1JJ2.out -infile_list  replace_scoreline_WITH_SHIFT_STATS_UNIQUE_H5_Feb_10_SWA.out  replace_scoreline_WITH_SHIFT_STATS_UNIQUE_H5_Dec_01_1JJ2.out   

			##################################################################################

		
			(shift_RMSD_col_name, shift_SCORE_col_name, num_shift_data_col_name) = get_shift_col_names(silent_file_type)

			keep_cols_list=["score", "O_loop_rmsd", shift_RMSD_col_name, shift_SCORE_col_name,  num_shift_data_col_name]

			ALL_PRE_silent_file_list=[]

			for data_src_ID in range(len(match_data_src_list)):

				data_src=match_data_src_list[data_src_ID]

				for ii in range(len(data_src["PRE_silent_file_list"])):

					PRE_silent_file=data_src["PRE_silent_file_list"][ii]

					submit_subprocess("cp %s %s" %(PRE_silent_file, basename(PRE_silent_file)))

					PRE_silent_file=basename(PRE_silent_file)

					################################################################################################################
					tag_prestring=""

					if(data_src["algorithm"]=="SWA_denovo"):
						tag_prestring="SWA_"
					elif(data_src["algorithm"]=="FARFAR_denovo"):
						tag_prestring="1JJ2_"
					else:
						error_exit_with_message("Invalid data_src[\"algorithm\"]=%s" %(data_src["algorithm"]))

					new_PRE_silent_file="replace_scoreline_%s" %(PRE_silent_file)

					replace_silent_scoreline_func(keep_cols_list, tag_prestring, PRE_silent_file, new_PRE_silent_file)

					PRE_silent_file=new_PRE_silent_file

					################################################################################################################
					if(silent_file_type=="PRE_MINIMIZE"):

						new_PRE_silent_file="%.1fX_NEW_CHEM_SHIFT_%s" %(int(pre_minimize_chem_shift_weight), PRE_silent_file)

						if(exists(new_PRE_silent_file)): error_exit_with_message("new_PRE_silent_file (%s) already exist!" %(new_PRE_silent_file))

						print 
						print "reweight PRE_MINIMIZE infile=(%s) | outfile=(%s) | pre_minimize_chem_shift_weight=(%s)" %(PRE_silent_file, new_PRE_silent_file, pre_minimize_chem_shift_weight)


						command ="reweight_chem_shift_score.py -infile %s -outfile %s " %(PRE_silent_file, new_PRE_silent_file)
						command+="-new_chem_shift_weight %s -use_rosetta_columns false " %(pre_minimize_chem_shift_weight)
						command+="> LOG_reweight_chem_shift_score_%s" %(PRE_silent_file)

						submit_subprocess(command)

						assert_is_valid_non_empty_silent_file(new_PRE_silent_file)

						PRE_silent_file=new_PRE_silent_file
					################################################################################################################

					data_src["PRE_silent_file_list"][ii]=PRE_silent_file
					ALL_PRE_silent_file_list.append(PRE_silent_file)
					################################################################################################################

			PRE_silent_file="PRE_silent_file.out"

			if(len(ALL_PRE_silent_file_list)==0): error_exit_with_message("len(ALL_PRE_silent_file_list)==0");
			if(exists(PRE_silent_file)): error_exit_with_message("PRE_silent_file (%s) already exist!" %(PRE_silent_file))
				
			concatenate_outfiles(infile_list=ALL_PRE_silent_file_list, outfile=PRE_silent_file, add_file_num_to_tag=False)

			assert_is_valid_non_empty_silent_file(PRE_silent_file)
			################################################################################################################
			recal_silent_file="calc_RMSD_to_native.out"

			recalculate_silent_file([PRE_silent_file], native_pdb, common_args_file, recal_silent_file)

			PRE_silent_file=recal_silent_file

			new_PRE_silent_file="replace_scoreline_%s" %(PRE_silent_file)

			spec_keep_cols_list=copy.deepcopy(keep_cols_list)
			spec_keep_cols_list.extend(["NEW_Full_L_rmsd-RMSD_V1_TO_NAT", "NEW_NAT_rmsd-RMSD_V2_TO_NAT"])
			keep_cols_list.extend(["RMSD_V1_TO_NAT", "RMSD_V2_TO_NAT"]) 

			replace_silent_scoreline_func(spec_keep_cols_list, "", PRE_silent_file, new_PRE_silent_file)

			PRE_silent_file=new_PRE_silent_file
		
			################################################################################################################
			best_SHIFT_RMSD_pdb=extract_best_model(PRE_silent_file, shift_RMSD_col_name)

			recal_silent_file="calc_RMSD_to_best_SHIFT_RMSD_MODEL.out"

			recalculate_silent_file([PRE_silent_file], best_SHIFT_RMSD_pdb, common_args_file, recal_silent_file)

			PRE_silent_file=recal_silent_file

			new_PRE_silent_file="replace_scoreline_%s" %(PRE_silent_file)

			spec_keep_cols_list=copy.deepcopy(keep_cols_list)
			spec_keep_cols_list.extend(["NEW_Full_L_rmsd-RMSD_V1_TO_BEST_SHIFT", "NEW_NAT_rmsd-RMSD_V2_TO_BEST_SHIFT"])
			keep_cols_list.extend(["RMSD_V1_TO_BEST_SHIFT", "RMSD_V2_TO_BEST_SHIFT"])  

			replace_silent_scoreline_func(spec_keep_cols_list, "", PRE_silent_file, new_PRE_silent_file)

			PRE_silent_file=new_PRE_silent_file

			################################################################################################################
			best_ENERGY_pdb=extract_best_model(PRE_silent_file, "score")

			recal_silent_file="calc_RMSD_to_best_ENERGY_MODEL.out"

			recalculate_silent_file([PRE_silent_file], best_ENERGY_pdb, common_args_file, recal_silent_file)

			PRE_silent_file=recal_silent_file

			new_PRE_silent_file="replace_scoreline_%s" %(PRE_silent_file)

			spec_keep_cols_list=copy.deepcopy(keep_cols_list)
			spec_keep_cols_list.extend(["NEW_Full_L_rmsd-RMSD_V1_TO_BEST_ENERGY", "NEW_NAT_rmsd-RMSD_V2_TO_BEST_ENEGY"])
			keep_cols_list.extend(["RMSD_V1_TO_BEST_ENERGY", "RMSD_V2_TO_BEST_ENEGY"])  

			replace_silent_scoreline_func(spec_keep_cols_list, "", PRE_silent_file, new_PRE_silent_file)

			PRE_silent_file=new_PRE_silent_file
			################################################################################################################

			os.chdir( working_foldername )

			if(exists(FINAL_CALC_silent_file)): error_exit_with_message("FINAL_CALC_silent_file (%s) already exist!")

			submit_subprocess("cp %s/%s %s" %(TEMP_FOLDER, PRE_silent_file, FINAL_CALC_silent_file))
	
			submit_subprocess("rm -r %s" %(TEMP_FOLDER))




print "------------------------------------------------------------------------------------------"
print "------------------------------------------------------------------------------------------"
print "GLOBAL: create_benchmark_table.py (script_algorithm=%s) RUN SUCCESSFULLY COMPLETED!" %(script_algorithm)
print "------------------------------------------------------------------------------------------"
print "------------------------------------------------------------------------------------------"






