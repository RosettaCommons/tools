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


from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
from SWA_dagman_python.parser.SWA_parse_benchmark import parse_benchmark_job_file
######################################################################



#create_benchmark_job_files.py -algorithm SWA -nstruct 100 -native_rmsd_screen True -job_list_file ../CS_benchmark_rna_list.txt  > outfile.txt 2> errfile.txt 
#create_benchmark_job_files.py -algorithm SWA -nstruct 5000 -native_rmsd_screen False -job_list_file ../CS_benchmark_rna_list.txt    > outfile.txt 2> errfile.txt 
#create_benchmark_job_files.py -algorithm SWA -nstruct 1000 -native_rmsd_screen False -job_list_file ../CS_benchmark_rna_list.txt    > outfile.txt 2> errfile.txt 


#create_benchmark_job_files.py -algorithm FARFAR -nstruct 250000 -nstruct_per_DAG 50000 -frac_struct_kept 0.20  -job_list_file  ../CS_benchmark_rna_list.txt    > outfile.txt 2> errfile.txt  
#create_benchmark_job_files.py -algorithm FARFAR -nstruct 250000 -nstruct_per_DAG 50000 -frac_struct_kept 0.50  -job_list_file  ../CS_benchmark_rna_list.txt    > outfile.txt 2> errfile.txt  
#create_benchmark_job_files.py -algorithm FARFAR -nstruct 250000 -nstruct_per_DAG 50000 -frac_struct_kept 0.02  -job_list_file  ../Oct_1_2011_benchmark_rna_list.txt  > outfile.txt 2> errfile.txt  

#~/minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/Nov_4_benchmark_rna_list.txt 
#~/minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/Oct_20_benchmark_rna_list.txt

######################################################################
main_folder=os.path.abspath("main_folder/")

if(exists(main_folder)): 
	print "%s already exist...removing..." %(main_folder)
	submit_subprocess("rm -r %s " %(main_folder) )
 
submit_subprocess( "mkdir %s" %(main_folder) )

os.chdir( main_folder)
######################################################################

job_list_file= parse_options( argv, "job_list_file", "../benchmark_rna_list.txt" )
nstruct= parse_options( argv, "nstruct", 1000 )
nstruct_per_DAG= parse_options( argv, "nstruct_per_DAG", 50000 ) #this option only effect FARFAR
frac_struct_kept= parse_options( argv, "frac_struct_kept", 0.02) #this option only effect FARFAR
native_rmsd_screen= parse_options( argv, "native_rmsd_screen", "False" )
num_slave_nodes=parse_options( argv, "num_slave_nodes", 500 ) #this is slave_jobs per particular motif.
algorithm=parse_options( argv, "algorithm", "SWA" ) 
star_mode= parse_options( argv, "star_mode", "True" )

dope_in_native_torsions=parse_options(argv, "dope_in_native_torsions", "False")
native_torsions_only=parse_options(argv, "native_torsions_only", "False")

if(dope_in_native_torsions==True and native_torsions_only==True): error_exit_with_message("dope_in_native_torsions and native_torsions_only cannot both be True!")

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s " %(list_to_string(argv) ) )

if(algorithm!="SWA" and algorithm!="FARFAR"): error_exit_with_message("Invalid algorithm (%s)! choices: (1) SWA (2) FARFAR" %(algorithm))  

if(algorithm=="FARFAR" and native_rmsd_screen==True): error_exit_with_message("native_rmsd_screen not supported in FARFAR" )  

native_rmsd_screen_str="False"
if(native_rmsd_screen): native_rmsd_screen_str="True"



job_list=parse_benchmark_job_file(job_list_file, setup_folders_and_files=True, star_mode=star_mode, verbose=True)


print "------------------------------------------------------------------------------------------"
print "------------------------------------------------------------------------------------------"
print "------------------------------------------------------------------------------------------"
print "------------------------------------------------------------------------------------------"


for job in job_list:

	print "create_dag_job_files for job " , job['folder_name']

	os.chdir( job['folder_name'])

	command=""
	if(algorithm=="SWA"):
		command= "create_dag_job_files.py -native_rmsd_screen %s "  %(native_rmsd_screen_str)
	else:
		command= "setup_FARFAR_job_files.py " 
		command+= "-nstruct_per_DAG %d " %(nstruct_per_DAG)
		command+= "-frac_struct_kept %f " %(frac_struct_kept)
		if(dope_in_native_torsions): command+= "-dope_in_native_torsions True"
		if(native_torsions_only):    command+= "-native_torsions_only True"


	if(job['native_pdb']=='fasta'): #new option Nov 22, 2010 pass in fasta instead of native_pdb
		command+=" -fasta_file %s " %(job['native_pdb'])
	else:
		command+=" -native_pdb %s " %(job['native_pdb'])

	if(job.has_key("nstruct")): error_exit_with_message("Please specify nstruct as a direct option to create_benchmark_job_files.py instead of to the job_file")

	command+=" -nstruct %d " %(nstruct)

	if(job['cutpoint_open']!="0"):
		command+=" -cutpoint_open %s "	%(list_to_string( job['cutpoint_open'].split("-") ) )

	if(job.has_key('num_slave_nodes')):
		command+= " -num_slave_nodes %d " %(int(job['num_slave_nodes'])) 
	else:
		command+= " -num_slave_nodes %d " %(num_slave_nodes) 

	if(job['bulge_res']!="0"): 
		command+=" -native_virtual_res %s " %(list_to_string( job['bulge_res'].split("-") ) )

	command+= " " + job['zextra'] + " "


	print command
	submit_subprocess(command)

	os.chdir( main_folder)

print "SUCCESSFULLY COMPLETED!"
