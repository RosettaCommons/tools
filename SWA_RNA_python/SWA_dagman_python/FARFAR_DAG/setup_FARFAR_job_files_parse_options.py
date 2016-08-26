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
from SWA_dagman_python.utility.torsion_database_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
from SWA_dagman_python.parser import SWA_parse_rosetta_options
######################################################################


#####################################################################################################################################################################################

def import_sequence(argv, verbose):

	input_sequence=parse_options( argv, "input_sequence", "" ).lower()
	fasta_file= parse_options( argv, "fasta_file", "" )
	native_pdb= parse_options( argv, "native_pdb", "" )

	if(native_pdb!="" and exists(native_pdb)==False): error_exit_with_message("native_pdb (%s) doesn't exist!!" %(native_pdb)) 

	if(fasta_file!="" and exists(fasta_file)==False): error_exit_with_message("fasta_file (%s) doesn't exist!!" %(fasta_file)) 

	num_input_sequence=0

	if(fasta_file!=""): num_input_sequence+=1
	if(native_pdb!=""): num_input_sequence+=1
	if(input_sequence!=""): num_input_sequence+=1

	if(num_input_sequence!=1): error_exit_with_message("num_input_sequence=%d!=1" %(num_input_sequence)) 

	#####create the fasta file if it doesn't exist##########
	if(fasta_file!=""):
		fasta_file=os.path.abspath(fasta_file)
	elif(native_pdb!=""):
		fasta_file=os.path.abspath("./fasta")
		submit_subprocess('SWA_pdb2fasta.py %s > %s' %(native_pdb,fasta_file))
	elif(input_sequence!=""):
		fasta_file=os.path.abspath("./fasta")
		FASTA_FILE=open( fasta_file, 'w')
		FASTA_FILE.write( ">input_sequence\n" )
		FASTA_FILE.write( "%s\n" %(input_sequence) )
		FASTA_FILE.close()
	else:
		error_exit_with_message("Cannot create fasta_file!")

	sequence = open( fasta_file ).readlines()[1][:-1] 
	if(verbose): print "sequence= %s " %(sequence)

	return (sequence, fasta_file, native_pdb)

#####################################################################################################################################################################################
def import_chunk_file_list(argv):

	chunk_file_list_method_1=parse_options( argv, "chunk_file_list", [""]) 

	chunk_file_list_method_2=parse_options( argv, "s", [""]) 

	if( chunk_file_list_method_1!=[""] and chunk_file_list_method_2!=[""] ):
		error_exit_with_message("User input chunk_files through both the -chunk_file_list and -s option")

	chunk_file_list=[""]

	if( chunk_file_list_method_1!=[""]): chunk_file_list=chunk_file_list_method_1

	if( chunk_file_list_method_2!=[""]): chunk_file_list=chunk_file_list_method_2

	print "imported_chunk_file_list= %s" %(chunk_file_list)

	return chunk_file_list


#####################################################################################################################################################################################

def import_nstruct_num_nodes(argv):
	###########Number structures to model and the number of CPU node available##########

	total_nstruct=parse_options(argv, "nstruct", 250000) 
	nstruct_per_DAG=parse_options(argv, "nstruct_per_DAG", 50000)
	frac_struct_kept=parse_options(argv, "frac_struct_kept", 0.02)
	num_slave_nodes=parse_options( argv, "num_slave_nodes", 500)

	if((total_nstruct % nstruct_per_DAG)!=0):  error_exit_with_message('total_nstruct (%d) not divisible by nstruct_per_DAG (%d)!' %(total_nstruct, nstruct_per_DAG ) )

	num_DAG=int(total_nstruct/nstruct_per_DAG)

	nstruct_per_node=int(nstruct_per_DAG/num_slave_nodes)  #this is per slave node
	if((nstruct_per_DAG%num_slave_nodes)!=0): nstruct_per_node+=1

	num_struct_kept=int((frac_struct_kept*nstruct_per_DAG)+0.1)

	print "------------------------------------------------------------------------------"
	print "total_nstruct=%d, nstruct_per_DAG=%d, frac_struct_kept=%f" %(total_nstruct, nstruct_per_DAG, frac_struct_kept)
	print "num_slave_nodes=%d, num_DAG= %d, nstruct_per_node=%d, num_struct_kept=%d" %(num_slave_nodes, num_DAG, nstruct_per_node, num_struct_kept) 
	print "------------------------------------------------------------------------------"

	return (num_slave_nodes, num_DAG, nstruct_per_node, num_struct_kept)


#####################################################################################################################################################################################


def import_torsion_database(argv, sample_res_list, native_pdb, sequence):

	torsion_database= parse_options( argv, "torsion_database", "1jj2.torsions") 

	excise_segment_torsion_DB = parse_options( argv , "excise_segment_torsion_DB", [""]) 
	ensure_excise_segment_seq_match = parse_options( argv , "ensure_excise_segment_seq_match", "True") 

	dope_in_native_torsions=parse_options(argv, "dope_in_native_torsions", "False")
	native_torsions_only=parse_options(argv, "native_torsions_only", "False")

	torsion_database="%s/chemical/rna/%s"  %(get_rosetta_database_folder(), torsion_database)
	print "INPUTTED torsion_database= %s " %(torsion_database)
	if(PATH_exists(torsion_database)==False): error_exit_with_message("torsion_database (%s) doesn't exist! " %(torsion_database) )

	######Since will be modify the torsion_database...safest to just create a local copy!
	if(torsion_database==basename(torsion_database)): error_exit_with_message("torsion_database==basename(torsion_database)")
	submit_subprocess("cp %s %s " %(torsion_database, basename(torsion_database) ) )  #create local_copy of torsion_databse
	torsion_database=basename(torsion_database)
	print "LOCAL_COPY torsion_database= %s " %(os.path.abspath(torsion_database))

	##########Create modifed torsion_database########
	if(excise_segment_torsion_DB!=[""]): #excise out native or homology torsions.
		torsion_database= create_excised_torsion_database(torsion_database, excise_segment_torsion_DB, ensure_excise_segment_seq_match, sequence, copy.deepcopy(sample_res_list) )

	if(dope_in_native_torsions or native_torsions_only):
		torsion_database= extract_native_torsions(dope_in_native_torsions, native_torsions_only, torsion_database, native_pdb, copy.deepcopy(sample_res_list) )

	if(dirname(torsion_database)!=""): error_exit_with_message("torsion_database_file (%s) is not in the main directory!" %(torsion_database) )

	return torsion_database

#####################################################################################################################################################################################

def get_rosetta_args_option(argv):

	common_args= SWA_parse_rosetta_options.get_general_rosetta_common_args_option(argv)

	rosetta_args= common_args

	####################################################################################################################
	lores_scorefxn=parse_options(argv, "lores_scorefxn", "") 

	sample_boundary_torsions=parse_options(argv, "sample_boundary_torsions", "True") #May 24, 2011 TO DO!!: Set to false to reproduce Rhiju's parameters.

	allow_consecutive_bulges=parse_options(argv, "allow_consecutive_bulges", "False")

	force_bulge_res=parse_options(argv, "force_bulge_res", [-1])

	fragment_insertion_cycles=parse_options(argv, "fragment_insertion_cycles", 10000)

	freeze_chunks_during_monte_carlo=parse_options(argv, "freeze_chunks_during_monte_carlo", "False")

	randomize_ribose_puckers=parse_options(argv, "randomize_ribose_puckers", "False") ##Oct 1st, 2011

	fragment_match_type=parse_options(argv, "fragment_match_type", "MATCH_YR")
	#Jan 05, 2012: Settle on MATCH_YR (compare MATCH_YR and MATCH_EXACT on SWA_loop_benchmark and gave similar results)

	####################################################################################################################

	if(lores_scorefxn!=""): rosetta_args += ' -lores_scorefxn %s ' %(lores_scorefxn) #April 9th, 2011

	if(sample_boundary_torsions==False): rosetta_args += ' -sample_boundary_torsions false '

	if(allow_consecutive_bulges==True): rosetta_args += ' -allow_consecutive_bulges true '

	rosetta_args += ' -cycles %d ' %(fragment_insertion_cycles) 

	if(len(force_bulge_res) > 0): rosetta_args += ' -virtual_res %s ' %(list_to_string(force_bulge_res) )  

	if(freeze_chunks_during_monte_carlo==True): rosetta_args += '-freeze_chunks_during_monte_carlo true '

	if(fragment_match_type!=""): 
		rosetta_args +="-fragment_match_type %s " %(fragment_match_type)

	if(randomize_ribose_puckers==True): rosetta_args += "-randomize_ribose_puckers true "

	####Extra option in SWA that have no meaning in FARFAR###############################################################
	BLAH_BLAH_BLAH=parse_options(argv, "OLLM_chain_closure_only", "", Verbose=False)
	BLAH_BLAH_BLAH=parse_options(argv, "allow_dinucleotide_at_loop_chain_closure", "", Verbose=False)
	BLAH_BLAH_BLAH=parse_options(argv, "enforce_path_base_pairs", [""] , Verbose=False)
	BLAH_BLAH_BLAH=parse_options(argv, "extra_anti_chi_rotamer", "false" , Verbose=False)
	BLAH_BLAH_BLAH=parse_options(argv, "extra_syn_chi_rotamer", "false" , Verbose=False)
	BLAH_BLAH_BLAH=parse_options(argv, "sample_virt_ribose_in_sep_DAG", "False" , Verbose=False)
	BLAH_BLAH_BLAH=parse_options(argv, "clusterer_keep_pose_in_memory", "true" , Verbose=False)
	BLAH_BLAH_BLAH=parse_options(argv, "clusterer_optimize_memory_usage", "false" , Verbose=False)
	BLAH_BLAH_BLAH=parse_options(argv, "old_SWA_idealize_helix", "False", Verbose=False) 
	#####################################################################################################################


	return (common_args,rosetta_args)

