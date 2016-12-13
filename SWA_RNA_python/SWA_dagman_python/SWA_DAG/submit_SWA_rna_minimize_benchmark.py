#!/usr/bin/env python

######################################################################
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.utility.DAGMAN_util import create_generic_README_SUB

from SWA_dagman_python.parser.SWA_parse_options import parse_options, get_option_name_args_safe, option_name_exist
from SWA_dagman_python.parser.SWA_parse_benchmark import *

######################################################################

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

#submit_SWA_rna_minimize_benchmark.py -num_slave_nodes_per_job 50 -star_mode False  -force_field_file stepwise/rna/rna_hires_07232011_with_intra_base_phosphate.wts  -remove_virtual_res_variant_during_minimize true > log_submit_swa_minimize_benchmark.out 2> log_submit_swa_minimize_benchmark.err

#-perform_minimize false 
#rna_4X_CS_03022011_w_intra_terms.wts

#submit_SWA_rna_minimize_benchmark.py -num_slave_nodes_per_job 250 -star_mode True  -force_field_file rna/rna_4X_CS_03022011_w_intra_terms.wts -remove_virtual_res_variant_during_minimize true > log_submit_swa_minimize_benchmark.out 2> log_submit_swa_minimize_benchmark.err

#

job_submission_file="/Users/sripakpa/minirosetta/Rosetta_rna_file/CS_PAPER/CS_benchmark_rna_list.txt"
data_src_file="/Users/sripakpa/minirosetta/Rosetta_rna_file/CS_PAPER/BENCHMARK_TABLE/DATA_SRC/March_05_2012_CS_benchmark_data_src.txt"
motif_file="/Users/sripakpa/minirosetta/Rosetta_rna_file/CS_PAPER/BENCHMARK_TABLE/DATA_SRC/CS_benchmark_motif_list.txt"

if(exists(job_submission_file)==False): error_exit_with_message("job_submission_file (%s) doesn't exist!!" %(job_submission_file))  
if(exists(data_src_file)==False): error_exit_with_message("data_src_file (%s) doesn't exist!!" %(data_src_file))  
if(exists(motif_file)==False): error_exit_with_message("motif_file (%s) doesn't exist!!" %(motif_file))  

print "job_submission_file=%s" %(job_submission_file)
print "data_src_file=%s" %(data_src_file)
print "motif_file=%s" %(motif_file)

######################################################
main_folder=os.path.abspath("main_folder/")
if(exists(main_folder)):
	"main_folder (%s) already exist! ... removing" %(main_folder)
	submit_subprocess("rm -r %s" %(main_folder))

submit_subprocess("mkdir %s" %(main_folder))

os.chdir( main_folder )

######################################################
submit_subprocess("cp %s %s" %(job_submission_file, 	basename(job_submission_file)))
submit_subprocess("cp %s %s" %(data_src_file, 				basename(data_src_file)))
submit_subprocess("cp %s %s" %(motif_file, 					basename(motif_file)))

job_submission_file=basename(job_submission_file)
data_src_file			=basename(data_src_file)
motif_file					=basename(motif_file)



######################################################
START_argv=copy.deepcopy(argv)

allow_missing_data= parse_options( argv, "allow_missing_data", "True" )
star_mode= parse_options( argv, "star_mode", "False" )
force_field_file=parse_options( argv, "force_field_file", "" )
num_slave_nodes_per_job=parse_options( argv, "num_slave_nodes_per_job", 100 )
remove_virtual_res_variant_during_minimize = parse_options( argv, "remove_virtual_res_variant_during_minimize", "false" )
perform_minimize = parse_options( argv, "perform_minimize", "true" )

if(force_field_file==""): error_exit_with_message("User need to specified force_field_file option!")

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )


DB = get_rosetta_database_folder()

if(PATH_exists('%s/scoring/weights/%s' %(DB, force_field_file) )==False): 
	error_exit_with_message("force_field_file (%s) doesn't exist!" %("%s/scoring/weights/%s" %(DB, force_field_file) ))

######################################################

motif_list=parse_motif_file(motif_file)

(data_src_list, algorithm_list)=parse_benchmark_data_source_file(data_src_file, star_mode, parse_silent_file=False)

job_submission_list=parse_benchmark_job_file(job_submission_file, setup_folders_and_files=False, star_mode=False, verbose=False) #ALWAYS set star_mode to false!


print "------------------------------------------------------------------------------------------"
print "------------------------------------------------------------------------------------------"
print "------------------------------------------------------------------------------------------"
print "------------------------------------------------------------------------------------------"

for motif in motif_list:

	#############################################################################################


	(job, found_job)=find_matching_job(job_submission_list, motif)

	if(found_job==False): error_exit_with_message("Cannot file job with job['folder_name']=%s" %(motif))

	#############################################################################################

	for algorithm in algorithm_list:

		os.chdir( main_folder )

		print "------------------------------------------------------------------------------------------"
		print "\nsubmitting run for motif=%s | algorithm=%s:" %(motif, algorithm)
		print 

		#############################################################################################
		found_data_src=False
		data_src={}

		for curr_data_src in data_src_list:

			if(curr_data_src['folder_name']==motif and curr_data_src['algorithm']==algorithm): 
				if(found_data_src): error_exit_with_message("Already found_data_src for motif=(%s) | algorithm=(%s)" %(motif , algorithm))
				data_src=copy.deepcopy(curr_data_src)
				found_data_src=True

		if(found_data_src==False): 
			print "\nMissing data_src for motif=(%s) | algorithm=(%s) ... Skipping data_src!\n" %(motif , algorithm)
			continue
		#############################################################################################
		
		submit_folder="%s/%s" %(algorithm, motif)

		if(exists(submit_folder)): error_exit_with_message("submit_folder (%s) already exist!" %(submit_folder))

		submit_subprocess("mkdir -p %s" %(submit_folder))			
		os.chdir( submit_folder )

		chemical_shift_file=get_option_name_args_safe(job['zextra'], "BMRB_chemical_shift_file", "")
		if(chemical_shift_file==""): error_exit_with_message("Unable to extract BMRB_chemical_shift_file option from job['zextra']!")

		chemical_shift_file=job['PRISTINE_motif_folder_name'] + '/' + chemical_shift_file

		native_pdb=job['PRISTINE_motif_folder_name'] + '/' + job['native_pdb']

		common_args_file= "%s/%s/%s" %(data_src['main_folder'], data_src['folder_name'], data_src['common_args'])

		PRE_minimize_silent_file= "%s/%s/%s" %(data_src['main_folder'], data_src['folder_name'], data_src['PRE_minimize_silent_file'])

		if(exists(chemical_shift_file)==False): error_exit_with_message("chemical_shift_file (%s) doesn't exist!" %(chemical_shift_file ))
		if(exists(native_pdb)==False): error_exit_with_message("native_pdb (%s) doesn't exist!" %(native_pdb))

		if(exists(common_args_file)==False): error_exit_with_message("common_args_file (%s) doesn't exist!" %(common_args_file))
		if(exists(PRE_minimize_silent_file)==False): error_exit_with_message("PRE_minimize_silent_file (%s) doesn't exist!" %( PRE_minimize_silent_file))

		submit_subprocess("cp %s %s" %(chemical_shift_file, basename(chemical_shift_file) ) )
		submit_subprocess("cp %s %s" %(native_pdb, basename(native_pdb) ) )
		submit_subprocess("cp %s %s" %(common_args_file, basename(common_args_file)))
		submit_subprocess("cp %s %s" %(PRE_minimize_silent_file, basename(PRE_minimize_silent_file)))

		chemical_shift_file=basename(chemical_shift_file)
		native_pdb=basename(native_pdb)
		common_args_file=basename(common_args_file)
		PRE_minimize_silent_file=basename(PRE_minimize_silent_file)

		####common_args_file implicitly need the VDW_rep_screener PDBs and fasta#####
		common_args_string=safe_readlines( common_args_file )[0]

		fasta_file=get_option_name_args_safe(common_args_string, "fasta", "")
		VDW_rep_screen_info=get_option_name_args_safe(common_args_string, "VDW_rep_screen_info", [""])

		submit_subprocess("SWA_pdb2fasta.py %s > %s" %(native_pdb, fasta_file))

		if(len(VDW_rep_screen_info) % 3 != 0): error_exit_with_message("len(VDW_rep_screen_info) % 3 != 0")

		for n in range(len(VDW_rep_screen_info)):

			if(n % 3 != 0): continue

			VDW_rep_screen_PDB="%s/%s/%s" %(data_src['main_folder'], data_src['folder_name'], VDW_rep_screen_info[n])

			if(exists(VDW_rep_screen_PDB)==False): error_exit_with_message("VDW_rep_screen_PDB (%s) doesn't exist!" %(VDW_rep_screen_PDB))

			submit_subprocess("cp %s %s" %(VDW_rep_screen_PDB, basename(VDW_rep_screen_PDB) ) )

		########################################################################
		create_generic_README_SUB(num_slave_nodes_per_job) 

		SWA_rna_minimize_command=  get_PYEXE("SWA_DAG/SWA_rna_minimize.py")

		SWA_rna_minimize_command+= " -silent_file %s " %(PRE_minimize_silent_file)

		SWA_rna_minimize_command+= " -force_field_file %s " %(force_field_file)

		SWA_rna_minimize_command+= " -chemical_shift_file %s " %(chemical_shift_file)

		SWA_rna_minimize_command+= " -chemical_shift_res_list DEFAULT "

		SWA_rna_minimize_command+= " -common_args %s " %(common_args_file)

		SWA_rna_minimize_command+= " -native_pdb %s " %(native_pdb)

		SWA_rna_minimize_command+= " -num_slave_nodes 0 " #This signals that wrapper function have already created README_SUB.py

		SWA_rna_minimize_command+= " -LOCATION BIOX2 "

		SWA_rna_minimize_command+= " -remove_virtual_res_variant_during_minimize %s " %(remove_virtual_res_variant_during_minimize)

		SWA_rna_minimize_command+= " -perform_minimize %s " %(perform_minimize)
		
		if( option_name_exist(job['zextra'], "chemical_shift_H5_prime_mode") ): 
			SWA_rna_minimize_command+= " -chemical_shift_H5_prime_mode %s " %(get_option_name_args_safe(job['zextra'], "chemical_shift_H5_prime_mode", ""))

		SWA_rna_minimize_command+= " > SWA_rna_minimize_LOG.txt "

		##~/SWA_RNA_python/SWA_dagman_python/SWA_DAG/SWA_rna_minimize.py -silent_file WITH_SHIFT_STATS_Nov_10_FARFAR.out  -force_field_file rna/rna_4X_CS_03022011_w_intra_terms.wts    -chemical_shift_file renumbered_AUTO_UCAC_shift.bmrb   -chemical_shift_res_list DEFAULT  -common_args common_args_region_FINAL.out -native_pdb G1C8_UCAC_1s72.pdb -num_slave_nodes 500    -LOCATION BIOX2 -remove_virtual_res_variant_during_minimize true > SWA_rna_minimize_LOG.txt 

		print "SWA_rna_minimize_command=%s" %(SWA_rna_minimize_command) 

		if(exists("README_SETUP.py")): error_exit_with_message("README_SETUP.py already exist!")

		README_SETUP = open( "README_SETUP.py", 'w')

		README_SETUP.write( '#!/usr/bin/env python\n' )
		README_SETUP.write( 'from os import system\n' )
		README_SETUP.write( 'import string\n\n' )
		README_SETUP.write( 'command=\"%s\"\n' %(SWA_rna_minimize_command	) )
		README_SETUP.write( 'system(command)\n' )
		README_SETUP.close()

		#submit_subprocess(SWA_rna_minimize_command)


print "------------------------------------------------------------------------------------------"
print "------------------------------------------------------------------------------------------"
print "GLOBAL: #submit_SWA_rna_minimize_benchmark.py RUN SUCCESSFULLY COMPLETED!" 
print "------------------------------------------------------------------------------------------"
print "START_argv=%s" %(START_argv)
print "------------------------------------------------------------------------------------------"
print "------------------------------------------------------------------------------------------"






