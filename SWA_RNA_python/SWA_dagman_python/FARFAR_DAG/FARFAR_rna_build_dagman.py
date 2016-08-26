#!/usr/bin/env python

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
import os
from os.path import basename, dirname, exists, expanduser
from time import sleep
from sets import Set
import copy
######################################################################

from FARFAR_dag_util import *
from FARFAR_rna_build_dagman_util import *
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.utility.DAGMAN_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################


#FARFAR_rna_build_dagman.py -num_DAG 5 

copy_argv=copy.deepcopy(argv)

final_rebuild_bulge = parse_options( argv, "final_rebuild_bulge", "False" )

BMRB_chemical_shift_file = parse_options( argv, "BMRB_chemical_shift_file", "" )
native_pdb = parse_options( argv, "native_pdb", "" )

filter_low_RMSD=(native_pdb!="")

rebuild_bulge_low_RMSD = parse_options( argv, "rebuild_bulge_low_RMSD", "False" )

if( (filter_low_RMSD==False) and (rebuild_bulge_low_RMSD==True) ): 
	print "Warning setting rebuild_bulge_low_RMSD to False since filter_low_RMSD==False"
	rebuild_bulge_low_RMSD=False


#########################################################################################################

num_DAG= parse_options( argv, "num_DAG", 0 )
num_struct_kept = parse_options( argv, "num_struct_kept", 0 )
num_slave_nodes= parse_options( argv, "num_slave_nodes", 0 )

if(num_struct_kept <= 0): error_exit_with_message('user need to pass in a positive integer for num_struct_kept')
if(num_DAG <= 0 ): error_exit_with_message('user need to pass in a positive integer for num_DAG')

if(num_slave_nodes <= 0): error_exit_with_message('user need to pass in a positive integer for num_slave_nodes')

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )

#########################################################################################################

fid_dag = open( "rna_build.dag", 'w' ) #Where the dag commands are ouputted.
fid_dag.write("DOT dag.dot\n")

filterer_job_tag_list=[] #doesn't include filterer jobs done is prior submission

if(exists("%s" %(DAG_FILE_FOLDER))): 
	print "WARNING DAG_FILE_FOLDER (%s) already exist! removing..." %(DAG_FILE_FOLDER)
	submit_subprocess("rm -r %s" %(DAG_FILE_FOLDER) )
submit_subprocess("mkdir %s" %(DAG_FILE_FOLDER) )

second_previous_filterer_job_tag=""
previous_filterer_job_tag=""

for DAG_ID in range(num_DAG):

	reducer_outfile="reducer_silent_file_DAG_ID_%d" %(DAG_ID)

	filtered_RMSD_file=get_filtered_RMSD_file(DAG_ID)
	filtered_energy_file=get_filtered_energy_file(DAG_ID)

	if(filter_low_RMSD):
		if(exists(filtered_RMSD_file)!=exists(filtered_energy_file)):
			error_exit_with_message('exists(filtered_RMSD_file)(%s)!=exists(filtered_energy_file)(%s)' %( exists( filtered_RMSD_file ),exists(filtered_energy_file) ) )

	if(exists(filtered_energy_file)): continue		#Jan 14, 2012: A little hacky since it is possible that filtered_energy_file exist but the filterer_job is not completely done

	##########################################################################

	DAG_FOLDER=get_DAG_ID_folder(DAG_ID)

	outfile_dir= "%s/$(Process)/" %(DAG_FOLDER)
	sampler_job_tag="DAG_ID_%d_SAMPLER" %(DAG_ID)

	empty_reducer_done_signal_file="%s/DAG_ID_%d_EMPTY_REDUCER_DONE.txt" %(DAG_FOLDER, DAG_ID)

	##########################################################################
	if(exists(empty_reducer_done_signal_file)==False):
	
		mapper_outfiles=outfile_dir + "silent_file.out" #WARNING THIS ASSUME SPECIFIC FILENAME IN FARFAR_rosetta_submit.py

		command_line = open(ROSETTA_COMMAND_FILE).readlines()[0]

		command_line = command_line[:-1].replace( '-out:file:silent ', '-out:file:silent '+ outfile_dir )

		if(is_release_mode()==False): #This is becuase I usually perform FARFAR setup locally on my Mac machine before rsync to cluster for job submssion!
			command_line = command_line.replace( 'macosgcc', 'linuxgcc')

		#command_line = command_line.replace( 'Users', 'home') #Comment out on Jan 01, 2012
		#command_line = command_line.replace( '/home/rhiju',HOMEDIR) #Comment out on Jan 01, 2012
		#command_line = command_line.replace( '/home/sripakpa',HOMEDIR) #Comment out on Jan 01, 2012

		params_file="params"
		if(exists(params_file)==False): error_exit_with_message("params_file (%s) doesn't exist!" %(params_file) )

		ROSETTA_DENOVO_EXE=command_line.split()[0]
		ROSETTA_DENOVO_arguments=list_to_string(command_line.split()[1:])

		sampler_job_filename = '%s/%s.condor' %(DAG_FILE_FOLDER, sampler_job_tag)

		fid_dag.write('SCRIPT PRE %s   %s %s %s\n' % (sampler_job_tag, PRE_PROCESS_SCRIPT , num_slave_nodes , sampler_job_filename) )

		fid_dag.write('\nJOB %s %s\n' % (sampler_job_tag, sampler_job_filename) )

		if(second_previous_filterer_job_tag!=""):
			fid_dag.write('PARENT %s CHILD %s\n' %(second_previous_filterer_job_tag, sampler_job_tag) )
		    
		make_dag_job_submit_file( sampler_job_filename, ROSETTA_DENOVO_EXE, ROSETTA_DENOVO_arguments, '', mapper_outfiles, '', '')

		fid_dag.write('SCRIPT POST %s %s -done_signal_file %s \n' %(sampler_job_tag, GENERIC_EMPTY_REDUCER_SCRIPT, empty_reducer_done_signal_file) )

	##########################################################################

	filterer_job_tag="DAG_ID_%d_FILTERER" %(DAG_ID)
	filterer_job_filename= '%s/%s.condor' %(DAG_FILE_FOLDER ,filterer_job_tag)
	filterer_arguments=" -DAG_ID %d -num_struct_kept %d -filter_low_RMSD %s " %(DAG_ID, num_struct_kept , filter_low_RMSD) 

	make_dag_job_submit_file( filterer_job_filename, FARFAR_FILTERER_EXE, filterer_arguments, '', '%s %s' %(filtered_RMSD_file,filtered_energy_file) , '', '')


	fid_dag.write('\nJOB %s %s\n' % (filterer_job_tag, filterer_job_filename) )

	if(exists(empty_reducer_done_signal_file)==False):
		fid_dag.write('PARENT %s CHILD %s\n' %(sampler_job_tag, filterer_job_tag) )

	filterer_job_tag_list.append(filterer_job_tag)

	second_previous_filterer_job_tag=previous_filterer_job_tag
	previous_filterer_job_tag=filterer_job_tag


########################################Final concatenation job#################################################
final_filtered_RMSD_file=   "FINAL_filtered_RMSD.out"
final_filtered_energy_file= "FINAL_filtered_energy.out"
FINAL_FILTER_tag="FINAL_FILTER" 

print "filterer_job_tag_list=%s " %(list_to_string(filterer_job_tag_list))
print "final_filtered_RMSD_file=%s " %(final_filtered_RMSD_file)
print "final_filtered_energy_file=%s " %(final_filtered_energy_file)

if(exists(final_filtered_energy_file)):

	if(filter_low_RMSD and (exists(final_filtered_RMSD_file)==False)): error_exit_with_message("final_filtered_energy_file  exist but final_filtered_RMSD_file does not exist!" )

	if(len(filterer_job_tag_list)>0): error_exit_with_message(" len(filterer_job_tag_list )=>0! but final_filtered_energy_file already exist!" )

else:

	FINAL_FILTER_filename= '%s/%s.condor' %(DAG_FILE_FOLDER, FINAL_FILTER_tag)

	FINAL_FILTER_arguments=" -num_DAG %s -filter_low_RMSD %s " %(num_DAG, filter_low_RMSD) 
	FINAL_FILTER_arguments+=" -final_filtered_RMSD_file %s " %(final_filtered_RMSD_file)
	FINAL_FILTER_arguments+=" -final_filtered_energy_file %s " %(final_filtered_energy_file)

	FINAL_FILTER_mapper_outfiles=final_filtered_energy_file
	if(filter_low_RMSD): FINAL_FILTER_mapper_outfiles+=' ' +final_filtered_RMSD_file

	make_dag_job_submit_file( FINAL_FILTER_filename, FARFAR_FINAL_FILTER_EXE, FINAL_FILTER_arguments, '', FINAL_FILTER_mapper_outfiles, '', '')

	fid_dag.write('\nJOB %s %s\n' % (FINAL_FILTER_tag, FINAL_FILTER_filename) )
	fid_dag.write('PARENT %s CHILD %s\n' %(string.join(filterer_job_tag_list), FINAL_FILTER_tag) )

########################################FINAL REBUILD BUILGE########################################################

common_args_split=(safe_readlines( COMMON_ARGS_FILE )[0][:-1]).split()
if(len(common_args_split)==0): error_exit_with_message("len(common_args_split)==0")

######First need to cluster the silent_file since it is not feasible to rebuild the bulge for every structures#####

FINAL_CLUSTER_ENERGY_job_tag="FINAL_CLUSTER_ENERGY"
FINAL_CLUSTER_RMSD_job_tag="FINAL_CLUSTER_RMSD"

clustered_energy_outfile=setup_FINAL_CLUSTER_job(fid_dag, FINAL_CLUSTER_ENERGY_job_tag, FINAL_FILTER_tag, final_filtered_energy_file)

if(rebuild_bulge_low_RMSD): 
	#Dec 09, 2011: If will not rebulge_bulge, then no need to cluster!
	#Dec 09, 2011: Also noticd that this clustering usually lead to code clash if use Richardson_RNA_09 Database!
	clustered_RMSD_outfile=setup_FINAL_CLUSTER_job(fid_dag, FINAL_CLUSTER_RMSD_job_tag, FINAL_FILTER_tag, final_filtered_RMSD_file)

####################################################################################################################

if(final_rebuild_bulge):

	if(BMRB_chemical_shift_file!=""):

		if(exists(BMRB_chemical_shift_file)==False): error_exit_with_message("BMRMB_chemical_shift_file (%s) doesn't exist!" %(BMRB_chemical_shift_file) )

		#-------------------------------------------------------------------------------------------------------#

		rebuild_bulge_job_tag_with_CS="FINAL_ENERGY_REBUILD_BULGE_WITH_CHEM_SHIFT"

		job_submitted=setup_final_rebuild_bulge_dag_job_file(fid_dag, clustered_energy_outfile, rebuild_bulge_job_tag_with_CS, native_pdb, \
																				 							COMMON_ARGS_FILE, "FARFAR", BMRB_chemical_shift_file)

		if( (exists(clustered_energy_outfile)==False) and job_submitted ): fid_dag.write('PARENT %s CHILD %s\n' %(FINAL_CLUSTER_ENERGY_job_tag, rebuild_bulge_job_tag_with_CS) ) 

		#-------------------------------------------------------------------------------------------------------#

		if(rebuild_bulge_low_RMSD): ###Added IF condition on Dec 09, 2011
			rebuild_bulge_job_tag_with_CS="FINAL_RMSD_REBUILD_BULGE_WITH_CHEM_SHIFT"

			job_submitted=setup_final_rebuild_bulge_dag_job_file(fid_dag, clustered_RMSD_outfile, rebuild_bulge_job_tag_with_CS, native_pdb, \
																				 							COMMON_ARGS_FILE, "FARFAR", BMRB_chemical_shift_file)

			if( (exists(clustered_RMSD_outfile)==False) and job_submitted ): fid_dag.write('PARENT %s CHILD %s\n' %(FINAL_CLUSTER_RMSD_job_tag, rebuild_bulge_job_tag_with_CS) ) 

	###################################################################################################################	
	rebuild_bulge_job_tag="FINAL_ENERGY_REBUILD_BULGE"

	job_submitted=setup_final_rebuild_bulge_dag_job_file(fid_dag, clustered_energy_outfile, rebuild_bulge_job_tag, native_pdb, COMMON_ARGS_FILE, "FARFAR", "")

	if( (exists(clustered_energy_outfile)==False) and job_submitted ): fid_dag.write('PARENT %s CHILD %s\n' %(FINAL_CLUSTER_ENERGY_job_tag, rebuild_bulge_job_tag) )

	#-------------------------------------------------------------------------------------------------------#

	if(rebuild_bulge_low_RMSD): ###Added IF condition on Dec 09, 2011
		rebuild_bulge_job_tag="FINAL_RMSD_REBUILD_BULGE"

		job_submitted=setup_final_rebuild_bulge_dag_job_file(fid_dag, clustered_RMSD_outfile, rebuild_bulge_job_tag, native_pdb, COMMON_ARGS_FILE, "FARFAR", "")

		if( (exists(clustered_RMSD_outfile)==False) and job_submitted ): fid_dag.write('PARENT %s CHILD %s\n' %(FINAL_CLUSTER_RMSD_job_tag, rebuild_bulge_job_tag) )

	###################################################################################################################

fid_dag.close()
#################################################################################################################

print "Algorithm %s successfully RAN! " %( list_to_string(copy_argv) )
print "----------------------------------------------------------------------------------------------------------------------------"

