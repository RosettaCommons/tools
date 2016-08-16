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
from setup_FARFAR_job_files_util import *
from setup_FARFAR_job_files_parse_options import *
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.utility.DAGMAN_util import setup_chemical_shift_args
from SWA_dagman_python.parser.SWA_parse_options import parse_options, parse_seq_num_list_option
######################################################################


#......((...((..((.....))....))...(()).....)).......
#AAAAAAUGAAAGGAAAAUUGGGUUAAAACCCCUUUAAUUCCCCAUUUUUUU
#12345678901234567890123456789012345678901234567890
#0        1         2         3         4    


#setup_FARFAR_job_files.py -input_sequence AAAAAAUGAAAGGAAAAUUGGGUUAAAACCCCUUUAAUUCCCCAUUUUUUU -cutpoint_open 20 35 -dot_bracket_structure "......((...((..((.....))....))...(()).....))......."


#setup_FARFAR_job_files.py -input_sequence GCGAAAGC  -dot_bracket_structure "((....))"


#setup_FARFAR_job_files.py -input_sequence GGCUUAAGCC -chunk_res_string_list 1-2 9-10 ! 4-7 -cutpoint_open 5

#######################################################Nov 08, 2011 FULL LENGTH##############################################################

#setup_FARFAR_job_files.py -native_pdb model1anr_RNA_A.pdb -num_slave_nodes 500  -chunk_file_list   HIV_TAR_lower_helix_EXCLUDE_A6_U24.pdb   Mock_upper_element_DAS_FARFAR_01.pdb -chunk_res_string_list 1-5 25-29 ! 11-22 -fixed_res 1-4 12-21 26-29  -fragment_match_type MATCH_YR -force_field_file rna_hires_07232011_with_intra_base_phosphate.wts > LOG_setup_FARFAR_job_files.out 2> LOG_setup_FARFAR_job_files.err


#######################################################Sept 06, 2011 MORE TETRALOOP STUFFF##############################################################

#setup_FARFAR_job_files.py -native_pdb Trimmed_2r8s_RNA.pdb -cutpoint_open 6 15 -num_slave_nodes 500 -chunk_file_list No_platform_Trimmed_2r8s_RNA.pdb  -chunk_res_string_list  1-10 13-23 -fixed_res 1-10 13-23 -rna_torsion_potential_folder FINAL_April_28_OLD_syn_chi/ -native_alignment_res 1 23 > LOG_setup_FARFAR_job_files.txt


#######################################################Aug 24, 2011 Tetraloop Revision Stuff############################################################

#setup_FARFAR_job_files.py -native_pdb U15_A16_Feb_10_2011_SWA_best_energy.pdb -cutpoint_open 6 15 -num_slave_nodes 250 -native_virtual_res 12 19 -chunk_file_list 2r8s_base_EX_U10_A18_U19.pdb CU_AG_BP.pdb -chunk_res_string_list  1-9 20-22 ! 14-17 -fixed_res 1 6 7 8 15 16 21 22 -rna_torsion_potential_folder FINAL_April_28_OLD_syn_chi/ > LOG_setup_FARFAR_job_files.txt

#setup_FARFAR_job_files.py -native_pdb Feb_10_2011_SWA_best_energy.pdb -cutpoint_open 6 15 -num_slave_nodes 500 -native_virtual_res 12 19 -chunk_file_list 2r8s_base_EX_U10_A18_U19.pdb  ccgg_build_helix_test.pdb -chunk_res_string_list  1-9 20-22 ! 14-17 -fixed_res 1 6 7 8 15 16 21 22 -rna_torsion_potential_folder FINAL_April_28_OLD_syn_chi/ > LOG_setup_FARFAR_job_files.txt


#setup_FARFAR_job_files.py -native_pdb Trimmed_2r8s_RNA.pdb -cutpoint_open 6 15 -num_slave_nodes 150 -native_virtual_res 20 -chunk_file_list 2r8s_base_element.pdb mutate_U2_C3_G4_A5_3BP_GU_invert_chain.pdb  -chunk_res_string_list  1-10 19-23 ! 13-18 -fixed_res 1-10 14-17 19-23 -rna_torsion_potential_folder FINAL_April_28_OLD_syn_chi/ > LOG_setup_FARFAR_job_files.txt

###########################################################################################################################################################

#setup_FARFAR_job_files.py  -fasta_file extended_purine_riboswitch_fasta   -chunk_res_string_list 1-15 71-91 ! 12-74  -chunk_file_list May_23_V3_cluster_200_FINAL_filtered_energy.out  remove_U12_A74_1y26_RNA_A.pdb -force_field_file 05232011_REP_only_rna_loop_hires.wts -lores_scorefxn 05232011_REP_only_rna_lores.wts -allow_bulge_mode False -freeze_chunks_during_monte_carlo True > LOG.txt

#setup_FARFAR_job_files.py -input_sequence GGCUUGGAAGCC -chunk_res_string_list 1 2 11 12 ! 4 5 8 9 -cutpoint_open 5

#setup_FARFAR_job_files.py -fasta_file 2r8s_fasta -chunk_res_string_list 1-120 126-144 149-159 -chunk_file_list start_chunk.pdb -nstruct 50000 -fragment_insertion_cycles 100000


########

HOMEDIR = expanduser('~') 

main_folder=os.path.abspath(".")

START_argv=copy.deepcopy(argv)

(COMMON_ARGS, ROSETTA_ARGS)=get_rosetta_args_option(argv)

verbose=parse_options( argv, "verbose", "True") 

final_rebuild_bulge=parse_options( argv, "final_rebuild_bulge", "True") 

##############################################################################
ignore_chemical_shift=parse_options( argv, "ignore_chemical_shift", "False") 

BMRB_chemical_shift_file=parse_options( argv, "BMRB_chemical_shift_file", "") 

if(ignore_chemical_shift): 
	BMRB_chemical_shift_file=""
	BLAHBLAHBLAH=parse_options( argv, "chemical_shift_H5_prime_mode", "" )
else:
	setup_chemical_shift_args(argv) 

##############################################################################

(sequence, fasta_file, native_pdb)=import_sequence(argv, verbose)

total_res=len(sequence)
all_res=range(1,total_res+1)

(num_slave_nodes, num_DAG, nstruct_per_node, num_struct_kept)=import_nstruct_num_nodes(argv)

chunk_res_list=[]
cutpoint_open_list=parse_options( argv, "cutpoint_open", [-1])

long_loop_mode=False

#####################Import chunk_res_list################################

if(len(chunk_res_list)==0): #import method #1 (from chunk_res_string_list)..this allow for pseudoknots and/or chunk_file containing nucleotides from two different branches

	chunk_res_string_combined=parse_options( argv, "chunk_res_string_list", [""] ) #example: 1 2 9 10 ! 4-7

	if(chunk_res_string_combined!=[""]):
		print "Import chunk_res_list method #1: chunk_res_string_list"
		chunk_res_string_combined=list_to_string(chunk_res_string_combined)
		for chunk_res_string in chunk_res_string_combined.split("!"):
			chunk_res_list.append(parse_segment_string_list(chunk_res_string.split()))

if(len(chunk_res_list)==0): #import method #2 (from dot_bracket_structure)..this doesn't allow for pseudoknots

	dot_bracket_structure=parse_options( argv, "dot_bracket_structure", "" ) 
	if(dot_bracket_structure!=""):
		print "Import chunk_res_list method #2: dot_bracket_structure"
		chunk_res_list=parse_dot_bracket_structure(dot_bracket_structure, total_res)

if(len(chunk_res_list)==0): #import method #3

	long_loop_mode=parse_options( argv, "long_loop_mode", "False") 
		
	if(long_loop_mode):
		print "Import chunk_res_list method #3: assuming long_loop_mode"
		sample_loop_res=parse_segment_string_list( parse_options(argv, "sample_res", [""]) )
		chunk_res_list.append( list(Set(all_res)-Set(sample_loop_res) ) )

if(len(chunk_res_list)==0): #import method #4
	
	helix_at_each_cutpoint_open=parse_options( argv, "helix_at_each_cutpoint_open", "True") 	

	if(helix_at_each_cutpoint_open):
		print "Import chunk_res_list method #4: assuming helix_at_each_cutpoint_open!"	
		helix_length_list=parse_options( argv, "helix_length_list", [-1])
		chunk_res_list=get_chunk_res_list_from_cutpoint_open(cutpoint_open_list, helix_length_list, total_res)		

print "NUM chunks= %s, chunk_list_res=" %(len(chunk_res_list)), chunk_res_list


######################Import the chunk pdbs/silent_files list######################

chunk_file_list=import_chunk_file_list(argv)	
fixed_res_list=parse_seq_num_list_option(argv, "fixed_res" )
VDW_rep_screen_info_list=parse_options(argv, "VDW_rep_screen_info", [""])
apply_VDW_rep_delete_matching_res=parse_options(argv, "apply_VDW_rep_delete_matching_res", "False")

apply_VDW_screener_to_helices=parse_options(argv, "apply_VDW_screener_to_helices", "True") #Switch this to True on March 07, 2012.
#if(apply_VDW_screener_to_helices): error_exit_with_message("apply_VDW_screener_to_helices==True, this option is not yet supported")

native_alignment_res=parse_options(argv, "native_alignment_res", [0])

if(VDW_rep_screen_info_list==[""]): VDW_rep_screen_info_list=[]

######################miscellaneous res_lists (not required)########################

allow_bulge_res_list=parse_options(argv, "allow_bulge_res", [-1])
allow_bulge_mode=parse_options(argv, "allow_bulge_mode", "True") #Nov 4, 2010...automatically virtualize bulge nucleotides


##################################################################################
virtual_phosphate_list=sorted( [1] +map( lambda x : x+1 , cutpoint_open_list ) ) #Determine this inside Rosetta instead?

(sample_res_list, sample_bb_list)=get_sample_regions(chunk_res_list, virtual_phosphate_list, total_res)
sample_segment_list=get_segment_string_list(sample_res_list, cutpoint_open_list)
rmsd_res_list=copy.deepcopy(sample_res_list) #might want to determine this internally inside Rosetta (to include sample_bb_list as well?)

allow_bulge_res_list=get_allow_bulge_res_list(allow_bulge_res_list, sample_res_list, cutpoint_open_list, total_res)

native_alignment_res=get_native_align_res_list(long_loop_mode, sample_res_list, sample_bb_list, total_res, native_alignment_res)

print "sample_segment_list=", sample_segment_list
print_seq_num_list("sample_res_list=", sample_res_list)
print_seq_num_list("sample_bb_list=", sample_bb_list)
print_seq_num_list("rmsd_res_list=", rmsd_res_list)
print_seq_num_list("allow_bulge_res_list=", allow_bulge_res_list)
print_seq_num_list("native_alignment_res=", native_alignment_res)
print_seq_num_list("virtual_phosphate_list=", virtual_phosphate_list)

torsion_database=import_torsion_database(argv, sample_res_list, native_pdb, sequence)

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )

chunk_info_list=get_chunk_info_list(chunk_res_list, chunk_file_list, sequence, sample_res_list, cutpoint_open_list )

fixed_res_list=get_fixed_res_list(long_loop_mode, fixed_res_list, sample_res_list, sample_bb_list, total_res)

(jump_point_list, additional_cutpoints_close)=get_jump_points_within_chunk_segments(chunk_info_list, fixed_res_list) #Rare!!

jump_point_list.extend( get_jump_points_between_chunk_segments(chunk_info_list, sample_res_list, sample_bb_list, fixed_res_list) )


###############create params file########################## 
if(exists("params")): error_exit_with_message('the file "params" already exists')

PARAMS=open( "params", 'w')

for n in range(len(jump_point_list)):
	###H H A  --> Hoogsteen-Hoogsteen Anti...anything here is find except for Anti Watson Crick!
	PARAMS.write( 'OBLIGATE   PAIR %d %d H H A\n'  %(jump_point_list[n][0], jump_point_list[n][1]) )

for segment in sample_segment_list: #This will be deprecate in favor of SAMPLE_BASE, SAMPLE_BB and SAMPLE_RES
	PARAMS.write( 'ALLOW_INSERT %s\n' %(segment.replace("-", " " ) ) )

for cutpoint_open in cutpoint_open_list:
	PARAMS.write( 'CUTPOINT_OPEN %d\n' %(cutpoint_open) )	

for chunk_info in chunk_info_list:
	#TO DO: distinguish START CHUNK!!!
	PARAMS.write( 'CHUNK %s %s\n' %( chunk_info["filename"] ,list_to_string(chunk_info["res_list"]) ) )	

if(len(additional_cutpoints_close)>0):
	PARAMS.write( 'ADDITIONAL_CUTPOINTS_CLOSE %s\n' %(list_to_string(additional_cutpoints_close) ) )

PARAMS.close()

###############create ROSETTA_COMMAND file#################
rosetta_command_file=get_FARFAR_ROSETTA_COMMAND_file()

if(exists(rosetta_command_file)): error_exit_with_message('the rosetta_command_file (%s) already exists' %(rosetta_command_file))

ROSETTA_COMMAND = open ( rosetta_command_file, 'w')
ROSETTA_COMMAND.write('%s  -database %s ' %(get_rosetta_EXE("rna_denovo"), get_rosetta_database_folder()) )

if(native_pdb!=""): ROSETTA_COMMAND.write('-native %s ' %(native_pdb))

ROSETTA_COMMAND.write('-fasta %s ' %(basename(fasta_file)) )
ROSETTA_COMMAND.write('-params_file params ' )
ROSETTA_COMMAND.write('-nstruct %d ' %(nstruct_per_node) )

if(len(fixed_res_list)>0): ROSETTA_COMMAND.write('-fixed_res %s ' %(list_to_string(fixed_res_list) ) ) #move this to params file

ROSETTA_COMMAND.write('-virtual_phosphate_list %s ' %(list_to_string(virtual_phosphate_list) ) ) #move this to inside Rosetta
ROSETTA_COMMAND.write('-rmsd_res %s ' %(list_to_string(rmsd_res_list) ) )
ROSETTA_COMMAND.write('-native_alignment_res %s ' %(list_to_string(native_alignment_res) ) )
ROSETTA_COMMAND.write('-vall_torsions %s ' %(torsion_database) )
ROSETTA_COMMAND.write('%s ' %(ROSETTA_ARGS) )
ROSETTA_COMMAND.write('-heat true -close_loops true -minimize_rna true -out:file:silent silent_file.out ' )


if(allow_bulge_mode):
	ROSETTA_COMMAND.write('-allow_bulge_mode true ' )
	ROSETTA_COMMAND.write('-allow_bulge_res_list %s ' %(list_to_string(allow_bulge_res_list) ) )

if(apply_VDW_screener_to_helices):
	for chunk_info in chunk_info_list:
		if(chunk_info["VDW_screener"]!=[]):
			VDW_rep_screen_info_list.extend(chunk_info["VDW_screener"])

if( len(VDW_rep_screen_info_list)>0 ):
	check_valid_VDW_rep_screen_info_list(VDW_rep_screen_info_list)

	if((len(VDW_rep_screen_info_list)%3)!=0): 
		print "VDW_rep_screen_info_list=", VDW_rep_screen_info_list
		error_exit_with_message("len(VDW_rep_screen_info_list)%3!=0") 

	if(long_loop_mode and len(VDW_rep_screen_info_list)!=3): error_exit_with_message("long_loop_mode and len(VDW_rep_screen_info_list)!=3") 

	if( (long_loop_mode==False) and (apply_VDW_screener_to_helices==False) ):
		error_exit_with_message("len(VDW_rep_screen_info_list)>0 but long_loop_mode==False and apply_VDW_screener_to_helices==False!")

	VDW_rep_screen_string=' -VDW_rep_screen_info %s ' %(list_to_string(VDW_rep_screen_info_list) ) 

	if(apply_VDW_rep_delete_matching_res): 
		if(long_loop_mode==False): error_exit_with_message("apply_VDW_rep_delete_matching_res but long_loop_mode==False!")
		if(len(sample_segment_list)==0): error_exit_with_message("len(sample_segment_list)==0")
		VDW_rep_screen_string+=' -VDW_rep_delete_matching_res %s ' %(list_to_string(sample_segment_list) )  
	else:
		VDW_rep_screen_string+=' -VDW_rep_delete_matching_res false '  

	ROSETTA_COMMAND.write(VDW_rep_screen_string)
	COMMON_ARGS+=VDW_rep_screen_string

ROSETTA_COMMAND.close()

##########create the common_args and cluster_args file (use to setup SWA job_params in final clustering and rebuild_bulge step)#########################
CLUSTER_ARGS=""

COMMON_ARGS+= ' -fasta %s ' %(basename(fasta_file)) 
COMMON_ARGS+= ' -rmsd_res %s ' %( list_to_string(rmsd_res_list))
COMMON_ARGS+= ' -native_alignment_res %s ' %( list_to_string(native_alignment_res))
COMMON_ARGS+= ' -global_sample_res_list %s ' %(list_to_string(sample_res_list) )

clusterer_alignment_res=[]
if(len(jump_point_list)==0):
	clusterer_alignment_res=[list_to_string(all_res,"-")[1:]] #Align over every nucleotides in the structure (optimal alignment)!
	CLUSTER_ARGS+= " -quick_alignment False "
else:
	clusterer_alignment_res= map(lambda jump_point : "%d-%d" %(jump_point[0],jump_point[1]) , jump_point_list )
	CLUSTER_ARGS+= " -quick_alignment True "


COMMON_ARGS+= ' -alignment_res %s ' %(list_to_string(clusterer_alignment_res) )

#Don't really need -cutpoint_open (currently only parsed by DAG_rebuild_bulge.py but actually not even used there if algorithm==FARFAR!!)
if(len(cutpoint_open_list)>0): COMMON_ARGS+= ' -cutpoint_open %s ' %(list_to_string(cutpoint_open_list) )

if(len(fixed_res_list)>0):
	COMMON_ARGS+= ' -fixed_res %s ' %(list_to_string(fixed_res_list))
else:
	COMMON_ARGS+= ' -minimize_res %s ' %(list_to_string(all_res))

COMMON_ARGS+= ' -simple_full_length_job_params true '

#########################################
cluster_args_file=get_FARFAR_cluster_args_file()

if(exists(cluster_args_file)): error_exit_with_message("cluster_args_file (%s) already exist!" %(cluster_args_file))

FARFAR_CLUSTER_ARGS=open( cluster_args_file, 'w')
FARFAR_CLUSTER_ARGS.write( CLUSTER_ARGS )
FARFAR_CLUSTER_ARGS.close()

#########################################
common_args_file=get_FARFAR_common_args_file()

if(exists(common_args_file)): error_exit_with_message("common_args_file (%s) already exist!" %(common_args_file))

FARFAR_COMMON_ARGS=open( common_args_file, 'w')
FARFAR_COMMON_ARGS.write( COMMON_ARGS )
FARFAR_COMMON_ARGS.close()

#########################################



################Don't need these for now!################################
#COMMON_ARGS+= ' -sample_res %d ' 
#COMMON_ARGS+= ' -cutpoint_closed %d ' 
#COMMON_ARGS+= ' -input_res %s '
#COMMON_ARGS+= ' -input_res2 %s ' 
#########################################################################

#############################create the DAG SUBMISSION FILE##############################

FARFAR_README_SETUP=open( "FARFAR_README_SETUP", 'w')
FARFAR_README_SETUP.write( "%s -num_slave_nodes %d -num_DAG %d "%(get_PYEXE("FARFAR_DAG/FARFAR_rna_build_dagman.py"), num_slave_nodes, num_DAG) )
FARFAR_README_SETUP.write( "-num_struct_kept %d -native_pdb %s "  %(num_struct_kept, native_pdb ) )

if(final_rebuild_bulge==True):
	FARFAR_README_SETUP.write( "-final_rebuild_bulge %s " %(final_rebuild_bulge) )
	if(BMRB_chemical_shift_file!=""):
		FARFAR_README_SETUP.write( "-BMRB_chemical_shift_file %s " %(BMRB_chemical_shift_file) )

FARFAR_README_SETUP.write( " > dagman_setup.txt \n" )

FARFAR_README_SETUP.close()

README_SUB = open( "FARFAR_README_SUB", 'w')
README_SUB.write( 'rm master_log.out\n' )
README_SUB.write( 'rm master_log.err\n' )
README_SUB.write( 'bsub -W 144:0 -o master_log.out -e master_log.err %s  -j %d rna_build.dag\n' %(get_PYEXE("dagman/DAG_continuous.py"), num_slave_nodes) )
README_SUB.close()

#######################################################################################################################

print "-------------------------------------------------------------------------------------------"
print "Successfully RAN: %s" %(list_to_string(START_argv))
print "-------------------------------------------------------------------------------------------"


