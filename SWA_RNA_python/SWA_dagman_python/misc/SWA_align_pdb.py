#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
import os
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options, parse_seq_num_list_option
######################################################################

#SWA_align_pdb.py -static_pdb first_BP_10_AU_BP.pdb -moving_pdb  FIRST_BP_4BP_5_UG_3_wobble.pdb  -static_align_res 1 -moving_align_res 1 -alignment_RMSD_CUTOFF 99.99

#SWA_align_pdb.py -static_pdb  nano_corner_2pn4.pdb -moving_pdb S_000000.pdb  -static_align_res 3-9 14-15 -moving_align_res 3-9 14-15 -alignment_RMSD_CUTOFF 99.99

#SWA_align_pdb.py -static_pdb  rosetta_Sigel_39mer_NMR_model.pdb -moving_pdb 1JJ2_S_47.pdb -static_align_res 2-5 8-12 -moving_align_res 2-5 8-12 -alignment_RMSD_CUTOFF 99.99

#SWA_align_pdb.py -static_pdb C10_FOX_NMR_WITH_WC_loop_ucac.pdb  -moving_pdb  SWA_S_38.pdb -static_align_res 3-8 -moving_align_res 3-8 -alignment_RMSD_CUTOFF 99.99

#SWA_align_pdb.py -static_pdb aligned_RNA_09_S_5626.pdb  -moving_pdb  SWA_S_5.pdb  -static_align_res 1 -moving_align_res 1 -alignment_RMSD_CUTOFF 99.99

#SWA_align_pdb.py -static_pdb DAS_FARFAR_01.pdb -moving_pdb  HIV_TAR_upper_helix.pdb -static_align_res 1 10 -moving_align_res 3  6 -alignment_RMSD_CUTOFF 0.001

#SWA_align_pdb.py -static_pdb LOOP_E_354d_RNA_A.pdb -moving_pdb S_000000.pdb -static_align_res 2-10 13-21 -moving_align_res 2-10 13-21 -alignment_RMSD_CUTOFF 99.0

#SWA_align_pdb.py -static_pdb C1G10_1040_1049_hexaloop_3huw.pdb -moving_pdb  S_1828.pdb -static_align_res 2-6 8-9 -moving_align_res 2-6 8-9 -alignment_RMSD_CUTOFF 99.0

#SWA_align_pdb.py -static_pdb short_GGAC_1mis.pdb  -moving_pdb  S_000000.pdb -static_align_res 2-5 8-11 -moving_align_res 2-5 8-11 -alignment_RMSD_CUTOFF 99.0

#SWA_align_pdb.py -static_pdb 1lnt_fix_non_natural_base.pdb  -moving_pdb  S_000000.pdb -static_align_res 2 3 4 5 6 7 10 11 12 13 14 15 -moving_align_res 2 3 4 5 6 7 10 11 12 13 14 15 -alignment_RMSD_CUTOFF 99.0

#SWA_align_pdb.py -static_pdb  duplex_1tut_RNA_A.pdb -moving_pdb S_000000.pdb -static_align_res 3-7 12-16 -moving_align_res 3-7 12-16 -alignment_RMSD_CUTOFF 5.0 

#SWA_align_pdb.py -static_pdb WITH_WC_BP_2koc.pdb -moving_pdb S_000000.pdb -static_align_res 3 4 6 7 8 -moving_align_res 3 4 6 7 8 -alignment_RMSD_CUTOFF 5.0 

#SWA_align_pdb.py -static_pdb DAS_FARFAR_03.pdb -moving_pdb  HIV_TAR_upper_helix.pdb -static_align_res 1  10 -moving_align_res 3  6 -alignment_RMSD_CUTOFF 0.001

#SWA_align_pdb.py -static_pdb 1lnt_fix_non_natural_base.pdb  -moving_pdb  S_000001.pdb -static_align_res 2 3 4 5 6 7 10 11 12 13 14 15 -moving_align_res 2 3 4 5 6 7 10 11 12 13 14 15 -alignment_RMSD_CUTOFF 99.0

#SWA_align_pdb.py -static_pdb lower_solution_hermann_duplex.pdb  -moving_pdb  LOWER_daslab_ts001_3.pdb -static_align_res 3 4 11 12 -moving_align_res 3 4 11 12 -alignment_RMSD_CUTOFF 1.0

#SWA_align_pdb.py -static_pdb hermann_phase2_puzzle_RNA_A.pdb  -moving_pdb lower_12_3_S_0.pdb  -static_align_res 1 2 3 14 15 16 -moving_align_res 1 2 3 12 13 14 -alignment_RMSD_CUTOFF 0.001
#SWA_align_pdb.py -static_pdb hermann_phase2_puzzle_RNA_A.pdb  -moving_pdb upper_12_3_S_0.pdb  -static_align_res 9 10 11 6 7 8 -moving_align_res 1 2 3 12 13 14 -alignment_RMSD_CUTOFF 0.001

#1-3 4-5 ,res 1 of static to res 3 of moving...res 4 of static to res 5 of moving.

#SWA_align_pdb.py -static_pdb S_000001_FARFAR_Nov_19_cluster_1A_hexaloop.pdb -moving_pdb mutate_A1-C2-C3-G4-G5-U6_3_BP_AU_helix.pdb -static_align_res 1 10 -moving_align_res 3 4 -alignment_RMSD_CUTOFF 0.001

#SWA_align_pdb.py -static_pdb inner_strand_2PN3_Square_RNA.pdb -moving_pdb /Users/sripakpa/minirosetta/Rosetta_rna_file/Square_RNA/create_VDW_rep_screener/A4U5_upper_VDW_rep_screener.pdb -alignment_res_pairs 6-5

#SWA_align_pdb.py -static_pdb inner_strand_2PN3_Square_RNA.pdb -moving_pdb A4U5_upper_VDW_rep_screener.pdb -alignment_res_pairs 6-5 
#
#2PN3_Square_RNA.pdb              2PN3_Square_RNA_fasta       

     
##############################################################

static_pdb =  parse_options( argv, "static_pdb", "" ) 
moving_pdb_list =	  parse_options( argv, "moving_pdb", [""] )
tag_name= parse_options( argv, "tag_name" , "" )

align_only_over_base_atoms=parse_options( argv, "align_only_over_base_atoms", "true")

if(tag_name!=""): 
	if(len(moving_pdb_list)!=1): error_exit_with_message("user passed tag_name (%s) but len(moving_pdb_list)!=1" %(tag_name) )


alignment_res_pairs =	  parse_options( argv, "alignment_res_pairs", [""] )
alignment_RMSD_CUTOFF =	  parse_options( argv, "alignment_RMSD_CUTOFF", 0.1 ) #ensure  perfect alignment

static_align_res =	parse_seq_num_list_option( argv, "static_align_res" ) 
moving_align_res =	parse_seq_num_list_option( argv, "moving_align_res" ) 

if(alignment_res_pairs!=[""]): 

	for n in range(len(alignment_res_pairs)):
		alignment_res_pair=alignment_res_pairs[n].split("-")
		if(len(alignment_res_pair)!=2): error_exit_with_message("len(alignment_res_pair)!=2!, alignment_res_pairs[%d]=%s" %(n, alignment_res_pairs[n]) )
elif(len(static_align_res)!=0):

	if(len(static_align_res)!=len(moving_align_res)):
		error_exit_with_message("len(static_align_res)!=len(moving_align_res), static_align_res=%s, moving_align_res=%s" %(static_align_res, moving_align_res))

	alignment_res_pairs=[]

	for ii in range(len(static_align_res)):
		alignment_res_pairs.append( "%d-%d" %( static_align_res[ii],moving_align_res[ii]) )

	if(len(alignment_res_pairs)==0): error_exit_with_message("?? len(alignment_res_pairs)==0")

else:
	error_exit_with_message("User need to pass in alignment_res_pairs or (static_align_res and moving_align_res)")



#if(len(alignment_res_pairs[0])==0 and len(alignment_res_pairs)==1): error_exit_with_message("User need to pass in alignment_res_pairs!")
#if(len(moving_pdb_list[0])==0 and len(moving_pdb_list)==1): error_exit_with_message("User need to pass in moving_pdb_list!")

if(moving_pdb_list==[""]): error_exit_with_message("User need to pass in moving_pdb_list!")



if(static_pdb==""): error_exit_with_message("User need to pass in static_pdb")

if(dirname(static_pdb)!=""):
	static_pdb=os.path.abspath(static_pdb)

if(exists(static_pdb)==False): error_exit_with_message("static_pdb (%s) doesn't exist!" %(static_pdb))


print "moving_pdb_list= " , moving_pdb_list
for n in range(len(moving_pdb_list)):
	moving_pdb=moving_pdb_list[n]
	if(dirname(moving_pdb)!=""): moving_pdb=os.path.abspath(moving_pdb)
	if(exists(moving_pdb)==False): error_exit_with_message("moving_pdb (%s) doesn't exist!" %(moving_pdb))
	
if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_args=%s " %(list_to_string(argv) ) )


##############################################################

if(use_new_src_code()):
	command=get_rosetta_EXE("swa_rna_util")
	command += ' -algorithm align_pdbs'
else:
	command=get_rosetta_EXE("parin_test")
	command += ' -algorithm alignment_test'


command += ' -database %s' %(get_rosetta_database_folder())
command += ' -native ' + static_pdb
command += ' -s ' + list_to_string(moving_pdb_list)
command += ' -alignment_res_pairs %s' %(list_to_string(alignment_res_pairs))
command += ' -alignment_RMSD_CUTOFF %s' %(alignment_RMSD_CUTOFF)
command += ' -align_only_over_base_atoms %s ' %(align_only_over_base_atoms)
if(tag_name!=""):
	command += ' -tag_name %s ' %(tag_name)
command += ' > SWA_align_pdb_output.txt '

print command 
submit_subprocess( command )




