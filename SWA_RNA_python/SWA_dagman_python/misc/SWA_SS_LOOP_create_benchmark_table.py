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

from SWA_SS_LOOP_create_benchmark_table_util import *

from SWA_dagman_python.misc.SWA_cat_outfiles import concatenate_outfiles

from SWA_dagman_python.utility.SWA_util import *

from SWA_dagman_python.parser.SWA_parse_internal_arguments import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
from SWA_dagman_python.parser.SWA_parse_benchmark import parse_benchmark_table_data_file, parse_motif_file
######################################################################


###SEPT 29, 2011. NOTE TO SELF SWA_filter_outfile_wrapper actually doesn't generate exactly "max_n_struct" structures. This is due to equivalency of the energy score are the cutoff (HAVE TO FIX THIS BEFORE the next time that I run the create_benchmark_table.py code!)

#create_benchmark_table.py -algorithm recal_rmsd -quick False -recal_shifts False -star_mode True  > log_recal_rmsd.out 2> log_recal_rmsd.err

#create_benchmark_table.py -algorithm create_table -allow_missing_data False -recreate_file True -star_mode True   > create_table.out 2> create_table.err



job_list_file="/Users/sripakpa/minirosetta/Rosetta_rna_file/CS_PAPER/BENCHMARK_TABLE/DATA_SRC/Nov_05_2011_CS_benchmark_table_data_src.txt"
motif_file="/Users/sripakpa/minirosetta/Rosetta_rna_file/CS_PAPER/BENCHMARK_TABLE/DATA_SRC/CS_benchmark_motif_list.txt"

#job_list_file= "Aug_28_benchmark_table_data_source.txt"
#if(job_list_file==""): error_exit_with_message('job_list_file==""')  
#job_list_file="/Users/sripakpa/minirosetta/Rosetta_rna_file/benchmark_list/benchmark_table_data_source/" + job_list_file
#motif_file=parse_options( argv, "motif_file", "LOOP_Benchmark_motif_list.txt")
#motif_file="/Users/sripakpa/minirosetta/Rosetta_rna_file/benchmark_list/benchmark_table_data_source/"+ motif_file 

if(exists(job_list_file)==False): error_exit_with_message("job_list_file (%s) doesn't exist!!" %(job_list_file))  
if(exists(motif_file)==False): error_exit_with_message("motif_file (%s) doesn't exist!!" %(motif_file))  

print "job_list_file=%s" %(job_list_file)
print "motif_file=%s" %(motif_file)

allow_missing_data= parse_options( argv, "allow_missing_data", "False" )
recreate_file= parse_options( argv, "recreate_file", "True" )
recal_shifts= parse_options( argv, "recal_shifts", "True" )
algorithm= parse_options( argv, "algorithm", "create_table" )
quick= parse_options( argv, "quick", "False" )
star_mode= parse_options( argv, "star_mode", "False" )


if(algorithm!="create_table" and algorithm!="recal_rmsd"):  error_exit_with_message("Invalid algorithm: %s" %(algorithm)) 


if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )

#####################################################

motif_list=parse_motif_file(motif_file)

######################################################


(job_list, algorithm_list)=parse_benchmark_table_data_file(job_list_file, algorithm, star_mode)


print "------------------------------------------------------------------------------------------"
print "------------------------------------------------------------------------------------------"
print "------------------------------------------------------------------------------------------"
print "------------------------------------------------------------------------------------------"

start_folder=os.path.abspath(".")

if(algorithm=="create_table"):

	table_filename="TABLE.txt"

	if(GET_OPTIMIZED_NATIVE_MODE): table_filename="OPTIMIZED_NATIVE_RMSD_%s_TABLE.txt" %(int(100*OPTIMIZED_NATIVE_RMSD_CUTOFF) )

	TABLE = open( table_filename, 'w')
	TABLE.write( '%s\n' %(job_list_file) )

	print_column_names(TABLE, algorithm_list)

	for motif in motif_list:

		if(DOUBLE_LINE_SPACING): TABLE.write( '\n' )
		TABLE.write( '%-34s' %(motif) )

		for algorithm in algorithm_list:

			found_data=False
			spec_job={}		

			for job in job_list:

				if(job['folder_name']==motif and job['algorithm']==algorithm): 
					spec_job=copy.deepcopy(job)
					found_data=True
					break

			############################################				
			if(found_data==False):
				print "motif=%s algorithm=%s doesn't exist in the job_list" %( motif, algorithm)
				for n in range(len(col_type_list)):	TABLE.write( '%-*s' %(COLUMNS_WIDTH,'') )

				continue
			############################################


			print "------------------------------------------------------------------------------------------"
			print "------------------------------------------------------------------------------------------"
			print "extracting data for motif=%s algorithm=%s :" %( motif, algorithm)
			print 


			full_foldername=""
			if(spec_job["main_folder"][-6:]=="$SELF$"):
				full_foldername="%s/" %( spec_job["main_folder"][:-6])
			else:
				full_foldername="%s/%s/" %( spec_job["main_folder"], spec_job["folder_name"])

			print "spec_job= ", spec_job

			print "------------------------------------------------------------------------------------------"
			print "------------------------------------------------------------------------------------------"

			if(exists(full_foldername)==False):
				if(allow_missing_data==False):
					error_exit_with_message("full_foldername (%s) doesn't exist!" %(full_foldername) )
				else:
					output_empty_string_filler(TABLE, len(col_type_list))
					print "full_foldername (%s) doesn't exist for motif=%s algorithm=%s!" %( full_foldername, motif, algorithm)
					continue						

			os.chdir( full_foldername )

			(spec_job['silent_file_list'], good_energy_silent_file_list)=get_recalculated_silent_file_list(spec_job, GET_OPTIMIZED_NATIVE_MODE)


			if(len(spec_job['silent_file_list'])==0):
				if(allow_missing_data==False):  
					error_exit_with_message("len(spec_job['silent_file_list'])==0, \n spec_job['silent_file']=%s \n full_foldername=%s " %(spec_job['silent_file'], full_foldername) )
				else:
					TABLE.write( '%-*s' %(COLUMNS_WIDTH, '') )
					print "len(spec_job['silent_file_list'])==0 for motif=%s algorithm=%s!" %( motif, algorithm )


			if(len(good_energy_silent_file_list)==0):
				if(allow_missing_data==False):  
					error_exit_with_message("len(good_energy_silent_file_list)==0, \n spec_job['silent_file']=%s \n full_foldername=%s " %(spec_job['silent_file'], full_foldername) )
				else:
					TABLE.write( '%-*s' %(COLUMNS_WIDTH, '') )
					print "len(good_energy_silent_file_list)==0 for motif=%s algorithm=%s!" %( motif, algorithm)


			if(exists(spec_job["native_pdb"])==False): error_exit_with_message("job[\"native_pdb\"] (%s) doesn't exist!: " %(spec_job["native_pdb"]) ) 

			#########Create the cluster_top_energy_silent_file######################################################
			cluster_silent_file=""
			SCORE_cluster_silent_file=""

			if( ( "best_RMSD_cluster" in col_type_list) or ("best_ENERGY_cluster" in col_type_list ) or OUTPUT_TOP_ENERGY_CLUSTERS_MODE):

				cluster_silent_file=create_cluster_silent_file(spec_job, good_energy_silent_file_list, recreate_file)

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
						output_top_energy_clusters_stats(TABLE, SCORE_cluster_silent_file, cluster_silent_file, spec_job, recreate_file)

					else:
						output_best_energy_and_best_RMSD_cluster_stats(TABLE, SCORE_cluster_silent_file, cluster_silent_file, spec_job, col_type_list, recreate_file)

				else: #col_type=="best_score" or col_type=="lowest_rmsd"

					output_best_score_or_best_rmsd(TABLE, spec_job, col_type, recreate_file)

			os.chdir( start_folder)
			TABLE.flush()

		TABLE.write('\n')

	TABLE.close()
else: #recal_rmsd or best_energy_clusters

	for job in job_list:
		print "------------------------------------------------------------------------------------------"
		print "------------------------------------------------------------------------------------------"
		if(algorithm=="recal_rmsd"):
			print "recalculating rmsd for job: ", 
		else:
			error_exit_with_message("Invalid algorithm (%s) " %(algorithm))

		print job['algorithm'] , " ", job["folder_name"]	, ":"
		print
		print job
		print "------------------------------------------------------------------------------------------"
		print "------------------------------------------------------------------------------------------"

		full_foldername=""
		if(job["main_folder"][-6:]=="$SELF$"):
			full_foldername="%s/" %( job["main_folder"][:-6])
		else:
			full_foldername="%s/%s/" %( job["main_folder"], job["folder_name"] )
		#print "full_foldername= %s" %(full_foldername)
		if(exists(full_foldername)==False): error_exit_with_message("full_foldername (%s) doesn't exist " %(full_foldername) )
		os.chdir( full_foldername )

		#print "job['common_args'][-6:]=", job["common_args"][-6:]

		if(exists(job["native_pdb"])==False): error_exit_with_message("job[\"native_pdb\"] (%s) doesn't exist!: " %(job["native_pdb"]) ) 

		common_args_file=job["common_args"]

		get_PRE_silent_file_list(job, quick)



		########################################################
		num_filter_struct=2000 
		#Change to 2000 on March 18
		#Change to 1000 on March 16..
		#Around March 7 change to 100...not sure what the value was before that..WARNING THIS IS NOT ENOUGH FOR FARFAR!

		PRE_rmsd_col_name="V_loop_rms" 
		#->V_loop_rmsd is more robust since virtual_res_list might change. Was "O_loop_rmsd" before March 17,  #Was "loop_rmsd" before March 8, 2011
		#->For loop benchmark, better to optimize_loop_rmsd since in some earlier run I specify native_alignment_res, but in rather run, this will taken off.
		if(job['algorithm'].count("FARFAR")>0): PRE_rmsd_col_name="Full_L_rmsd"
		########################################################

		quick_individual_regions=False
		if(job['PRE_silent_file']=="INDIVIDUAL_REGIONS" and quick==True): quick_individual_regions=True 	

		all_filtered_RMSD_silent_file_list=[]
		all_filtered_score_silent_file_list=[]

		########################################################
		working_foldername="%s/" %(algorithm.upper())

		if(quick):  working_foldername="FILTERED_" + working_foldername

		suffix= "recal_"


		######################################################## 
		if(CHEMICAL_SHIFT_MODE): 

			#In CHEMICAL_SHIFT_MODE, the new PRE_silent_file with the chemical_shift stats is in the BENCHMARK_TABLE_CALC_CHEM_SHIFT/ folder!

			job["CHEM_shift_silent_file_list"]=[]

			if(recal_shifts):

				for PRE_silent_file in job["PRE_silent_file_list"]:

					calc_chem_shift_folder="BENCHMARK_TABLE_CALC_CHEM_SHIFT/"

					if(dirname(PRE_silent_file)!=""): calc_chem_shift_folder="%s/%s" %(dirname(PRE_silent_file), calc_chem_shift_folder)

					if(exists(calc_chem_shift_folder)):
						print "calc_chem_shift_folder (%s) already exist ... removing!" %(calc_chem_shift_folder)
						submit_subprocess( "rm -r %s " %(calc_chem_shift_folder) ) #Might as well delete the whole folder since will recal_shift


			for PRE_silent_file in job["PRE_silent_file_list"]:

				calc_chem_shift_folder="BENCHMARK_TABLE_CALC_CHEM_SHIFT/"

				if(dirname(PRE_silent_file)!=""): calc_chem_shift_folder="%s/%s" %(dirname(PRE_silent_file), calc_chem_shift_folder)

				submit_subprocess( "mkdir -p %s " %(calc_chem_shift_folder) )

				job["CHEM_shift_silent_file_list"].append( create_silent_file_WITH_chemical_shift_stats(job, PRE_silent_file, calc_chem_shift_folder, recal_shifts) )

		######################################################## 

		if(recreate_file): #This clean up working_foldername from prior run if recreate_file == True..BUT not the Chemical_shift_folder!

			for PRE_silent_file in job["PRE_silent_file_list"]: 

				specific_working_foldername=working_foldername

				if(dirname(PRE_silent_file)!=""): specific_working_foldername="%s/%s" %(dirname(PRE_silent_file), working_foldername)

				if(exists(specific_working_foldername)):
					print "specific_working_foldername (%s) already exist ... removing!" %(specific_working_foldername)
					submit_subprocess( "rm -r %s " %(specific_working_foldername) ) #Might as well delete the whole folder since will recreate the file

		########################################################

		for n in range(len(job["PRE_silent_file_list"])):

			PRE_silent_file=job["PRE_silent_file_list"][n]

			if(check_PRE_silent_file(PRE_silent_file)): continue

			specific_working_foldername=working_foldername

			if(dirname(PRE_silent_file)!=""): specific_working_foldername="%s/%s" %(dirname(PRE_silent_file), working_foldername)

			submit_subprocess( "mkdir -p %s " %(specific_working_foldername) )

			####Hacky..change to WITH_chemical_shift_silent after creating the working_folder!#####
			if(CHEMICAL_SHIFT_MODE): PRE_silent_file=job["CHEM_shift_silent_file_list"][n]

			if(quick==False): 

				recal_silent_file=specific_working_foldername + suffix + "full_"  + basename(PRE_silent_file)

				#recal_silent_file=append_to_basename(suffix + "full_" , PRE_silent_file)

				if(check_for_existing_recalculate_silent_file(recreate_file, recal_silent_file)): continue
		
				recalculate_silent_file([PRE_silent_file], job, common_args_file, recal_silent_file)

			else: 
				
				####################SCORE########################################
				#recal_filter_silent_file=append_to_basename(suffix + "filtered_score_" , PRE_silent_file) //Nov 06, 2011
				recal_filter_silent_file=specific_working_foldername + suffix + "filter_score_" + basename(PRE_silent_file)
				

				if(quick_individual_regions or check_for_existing_recalculate_silent_file(recreate_file, recal_filter_silent_file)==False): 

					#filter_outfile=append_to_basename("%s%s_max_n_struct_%d_" %(working_foldername, "score", num_filter_struct), PRE_silent_file ) //Nov 06, 2011 
					filter_outfile= "%s%s_max_n_struct_%d_%s" %(specific_working_foldername, "score", num_filter_struct, basename(PRE_silent_file ) )

					filter_command="SWA_filter_outfile_wrapper.py -infile %s -scorecol_name %s -max_n_struct %d -filter_outfile %s -remove_SCORE_file True" %(PRE_silent_file, "score", num_filter_struct, filter_outfile) 

					submit_subprocess(filter_command)	

					all_filtered_score_silent_file_list.append(os.path.abspath(filter_outfile))

					if(quick_individual_regions):
						all_filtered_score_silent_file_list.append(os.path.abspath(filter_outfile))
					else:
						if(check_for_existing_recalculate_silent_file(recreate_file, recal_filter_silent_file)==False): 
							recalculate_silent_file([filter_outfile], job, common_args_file, recal_filter_silent_file)

				#####################RMSD########################################
				#recal_filter_silent_file=append_to_basename(suffix + "filtered_RMSD_",  PRE_silent_file) //Nov 06, 2011
				recal_filter_silent_file=specific_working_foldername + suffix + "filter_RMSD_" + basename(PRE_silent_file)

				if(quick_individual_regions or check_for_existing_recalculate_silent_file(recreate_file, recal_filter_silent_file)==False): 

					#filter_outfile=append_to_basename("%s%s_max_n_struct_%d_" %(working_foldername, PRE_rmsd_col_name, num_filter_struct), PRE_silent_file ) //Nov 06, 2011
					filter_outfile= "%s%s_max_n_struct_%d_%s" %(specific_working_foldername, PRE_rmsd_col_name, num_filter_struct, basename(PRE_silent_file ) )

					filter_command="SWA_filter_outfile_wrapper.py -infile %s -scorecol_name %s -max_n_struct %d -filter_outfile %s -remove_SCORE_file True" %(PRE_silent_file, PRE_rmsd_col_name, num_filter_struct, filter_outfile) 

					submit_subprocess(filter_command)	

					if(quick_individual_regions):
						all_filtered_RMSD_silent_file_list.append(os.path.abspath(filter_outfile))
					else:
						if(check_for_existing_recalculate_silent_file(recreate_file, recal_filter_silent_file)==False): 
							recalculate_silent_file([filter_outfile], job, common_args_file, recal_filter_silent_file)


		########################################################################################################
		if(quick_individual_regions):

			if(len(all_filtered_RMSD_silent_file_list)==0): error_exit_with_message("len(all_filtered_RMSD_silent_file_list)==0!")
			if(len(all_filtered_score_silent_file_list)==0): error_exit_with_message("len(all_filtered_score_silent_file_list)==0!")


			#1. Combine all the filtered silent_files for all best RMSD and for all best ENERGY before recalculate silent_file
			#2. Important to use the same num_filter_struct as in filtering the individual files as to NOT MISS any data points that would exist in the non-quick mode.

			working_foldername="INDIVIDUAL_REGIONS/" + working_foldername

			if(recreate_file and exists(working_foldername)):
				print "working_foldername (%s) already exist ... removing!" %(working_foldername)
				submit_subprocess("rm -r %s " %(working_foldername) ) #Might as well delete the whole folder since will recreate the file

			submit_subprocess( "mkdir -p %s " %(working_foldername ) )

			##############################SCORE########################
			recalculate_score_silent_file=working_foldername +  "/recal_filtered_score_individual_regions.out"

			if( check_for_existing_recalculate_silent_file(recreate_file, recalculate_score_silent_file)==False ):

				print "all_filtered_score_silent_file_list=", all_filtered_score_silent_file_list
				cat_filtered_score_file="%s/%s" %(working_foldername, "cat_filtered_score_individual_regions.out")
				print "cat_filtered_score_file=", cat_filtered_score_file

				if(exists(cat_filtered_score_file)):
					print "cat_filtered_score_file (%s) already exist! ...removing..." %(cat_filtered_score_file)
					submit_subprocess("rm %s " %(cat_filtered_score_file))

				concatenate_outfiles(infile_list=copy.deepcopy(all_filtered_score_silent_file_list), outfile=cat_filtered_score_file) 				


				final_filter_score_outfile=append_to_basename("score_max_n_struct_%d_" %(num_filter_struct), cat_filtered_score_file ) 
				filter_command="SWA_filter_outfile_wrapper.py -infile %s -scorecol_name %s -max_n_struct %d -filter_outfile %s -remove_SCORE_file True" %(cat_filtered_score_file , "score", num_filter_struct, final_filter_score_outfile) 
				submit_subprocess(filter_command)	

				recalculate_silent_file([final_filter_score_outfile], job, common_args_file, recalculate_score_silent_file)

			##############################RMSD########################
			recalculate_RMSD_silent_file=working_foldername + "/recal_filtered_RMSD_individual_regions.out"

			if( check_for_existing_recalculate_silent_file(recreate_file, recalculate_RMSD_silent_file)==False ):

				print "all_filtered_RMSD_silent_file_list=", all_filtered_RMSD_silent_file_list
				cat_filtered_RMSD_file="%s/%s" %(working_foldername, "cat_filtered_RMSD_individual_regions.out")
				print "cat_filtered_RMSD_file=", cat_filtered_RMSD_file

				if(exists(cat_filtered_RMSD_file)):
					print "cat_filtered_RMSD_file (%s) already exist! ...removing..." %(cat_filtered_RMSD_file)
					submit_subprocess("rm %s " %(cat_filtered_RMSD_file))

				concatenate_outfiles(infile_list=copy.deepcopy(all_filtered_RMSD_silent_file_list), outfile=cat_filtered_RMSD_file) 				

				final_filter_RMSD_outfile=append_to_basename("%s_max_n_struct_%d_" %(PRE_rmsd_col_name, num_filter_struct), cat_filtered_RMSD_file ) 
				filter_command="SWA_filter_outfile_wrapper.py -infile %s -scorecol_name %s -max_n_struct %d -filter_outfile %s -remove_SCORE_file True" %(cat_filtered_RMSD_file , PRE_rmsd_col_name, num_filter_struct, final_filter_RMSD_outfile) 

				submit_subprocess(filter_command)	

				recalculate_silent_file([final_filter_RMSD_outfile], job, common_args_file, recalculate_RMSD_silent_file)

			##########################################################

		########################################################################################################


	os.chdir( start_folder)



print "------------------------------------------------------------------------------------------"
print "------------------------------------------------------------------------------------------"
print "GLOBAL: create_benchmark_table.py (algorithm=%s) RUN SUCCESSFULLY COMPLETED!" %(algorithm)
print "------------------------------------------------------------------------------------------"
print "------------------------------------------------------------------------------------------"






