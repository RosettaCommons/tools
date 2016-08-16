#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
import os
import copy

######################################################################
from SWA_util import *
from RNA_BP_and_BS_util import *
######################################################################

def convert_MC_to_LW_base_edge(edge, include_bifurcation_BP):

	####Refer to "RNA canonical and non-canonical base pairing types: a recognition method and complete repertoire" by the Major Group

	if(edge=="Ww" or edge=="Ws" or edge=="Wh"): return "W"

	if(edge=="Ss" or edge=="Sw"): return "S"

	if(edge=="Hh" or edge=="Hw"): return "H"

	if(edge=="C8"): return "N" #"H" #From 2002 paper: "We also introduced a special face, C8, for the C8-H8 donor group of the purines." NOT a polar Hbond!

	if(edge=="O2'"): return "N" #N to none

	if(edge=="O1P"): return "N" #N to none

	if(edge=="O2P"): return "N" #N to none

	if(include_bifurcation_BP): 
		#error_exit_with_message("include_bifurcation_BP mode not implemented yet!")
		if(edge=="Bs"): return "Bs"
		if(edge=="Bh"): return "Bh"
	else:
		if(edge=="Bs"): return "N" #N to none
		if(edge=="Bh"): return "N" #N to none

	#consistency check
	if(edge[0]=="W"): error_exit_with_message("A unaccounted Watson-crick edge!")
	if(edge[0]=="S"): error_exit_with_message("A unaccounted Sugar edge!")
	if(edge[0]=="H"): error_exit_with_message("A unaccounted Hoogsten edge!")
	if(edge[0]=="B"): error_exit_with_message("A unaccounted Bifurcation edge!")

	error_exit_with_message("Invalid base_pair edge (%s) " %(edge) )


	####OK deal with bifurcated Base-pair (Bh and Bs) later!


def extract_MC_annotate_data_func(rosetta_pdb_name, sample_segment, double_count=True, include_bifurcation_BP=False, Verbose=True):

	MC_Annotate="/Users/sripakpa/Program_labtop/MC_ANNOTATE/MC-Annotate"

	if( exists( rosetta_pdb_name )==False):  error_exit_with_message("rosetta_pdb_name (%s) doesn't exist!" %(rosetta_pdb_name) )

	rosetta_pdb_name=os.path.abspath(rosetta_pdb_name)

	###############################################################################################
	fasta_file = "fasta"

	submit_subprocess("SWA_pdb2fasta.py %s > %s" %(rosetta_pdb_name, fasta_file) )

	sequence = open( fasta_file  ).readlines()[1][:-1]

	SEQUENCE=sequence.upper()

	###############################################################################################

	if(sample_segment=="All"):
		start_res=1
		end_res=len(SEQUENCE)

	else:
		sample_segment=sample_segment.split('-')
		if(len(sample_segment)!=2): error_exit_with_message("len(sample_segment)!=2" )

		start_res=int(sample_segment[0])
		end_res=int(sample_segment[1])

		###consistency check########

		if(start_res<1):  error_exit_with_message("start_res is not a positive integer, start_res=%d" %(start_res) )
		if(end_res<1):  error_exit_with_message("end_res is not positive integer, end_res=%d" %(end_res) )

		if(start_res>=end_res): error_exit_with_message("start_res(%d) > end_res(%d)" %(start_res, end_res))

	###############################################################################################


	sample_res=''
	loop_res_list=[]

	for seq_num in range(start_res, end_res+1):
		sample_res+='%d ' %(seq_num)

		loop_res_list.append(seq_num)
		
	###############################################################################################

	foldername="extract_MC_ANNOTATE_data_%s_%4s_%4s" %( basename(rosetta_pdb_name).replace(".pdb",""), str(start_res).zfill(4), str(end_res).zfill(4))

	if( exists(foldername) ):   
		if(Verbose): print "Warning foldername (%s) already exist! ...removing!..." %(foldername) 
		submit_subprocess("rm -r %s " %(foldername) )


	submit_subprocess( 'mkdir %s' %(foldername)  )

	os.chdir( foldername)


	####################Convert Rosetta format PDB to standard format PDB##########################

	standard_pdb_name="STANDARD_%s" %(basename(rosetta_pdb_name))

	submit_subprocess("Rosetta_to_standard_PDB.py -s %s -output_pdb %s > LOG_Rosetta_to_standard_PDB.txt" %(rosetta_pdb_name, standard_pdb_name) )


	###############################################################################################
	MC_data_file="MC_data_" + standard_pdb_name[0:-4] + ".txt"

	MC_command="%s %s > %s" %(MC_Annotate , standard_pdb_name, MC_data_file)
	if(Verbose): print "MC_command=%s " %(MC_command)
	submit_subprocess( MC_command ) #Choose allow sending of multiple pdb?


	###############################################################################################

	if(exists(MC_data_file)==False): error_exit_with_message("MC_data_file (%s) doesn't exist!" %(MC_data_file) )

	MC_lines=open(MC_data_file).readlines()


	if(MC_lines[0]!="Residue conformations -------------------------------------------\n"): error_exit_with_message("MC_lines[0]!=\"Residue conformations -------------------------------------------\"")

	line_num=0

	while(True):

		line_num+=1

		line=MC_lines[line_num]

		line_list=line.split()

		if(line=="Adjacent stackings ----------------------------------------------\n"): break

		if(line.count("--")>0): error_exit_with_message("line.count(\"-\")>0 for line=%s" %(line) )

		seq_num=int(line_list[0][1:])

		MC_SEQ=line_list[2]

		if(seq_num!=line_num): error_exit_with_message("seq_num=(%s)!=(%s)=line_num" %(seq_num,line_num) )


		if(MC_SEQ!=SEQUENCE[seq_num-1]): error_exit_with_message("MC_SEQ=(%s)!=(%s)=SEQUENCE[seq_num-1]" %(MC_SEQ,SEQUENCE[seq_num-1] ) )

	Base_pair_list=[]
	Base_stack_list=[]


	while(True):

		line_num+=1

		line=MC_lines[line_num]

		line_list=line.split()

		if(line=="Non-Adjacent stackings ------------------------------------------\n"): continue

		if(line[:6]=="Number"): continue

		if(line=="Base-pairs ------------------------------------------------------\n"): break

		if(line.count("--")>0): error_exit_with_message("line.count(\"-\")>0 for line=%s" %(line) )

		seq_num_1=int(line_list[0].split("-")[0][1:])
		seq_num_2=int(line_list[0].split("-")[1][1:])

		if( (seq_num_1 not in loop_res_list) and (seq_num_2 not in loop_res_list) ): continue # not an interaction involving the loop!

		if(seq_num_1==seq_num_2): error_exit_with_message("seq_num_1=(%s)==(%s)=seq_num_2!", (seq_num_1,seq_num_2) )

		interaction=list_to_string(line_list[2:], first_separator=False)

		face1=""
		face2=""

		if(interaction.count("upward")>0): 
			face1="3"
			face2="5"
		elif	(interaction.count("downward")>0): 
			face1="5"
			face2="3"
		elif	(interaction.count("inward")>0): 
			face1="3"
			face2="3"
		elif	(interaction.count("outward")>0): 
			face1="5"
			face2="5"
		else:
			error_exit_with_message("Invalid base-stack interaction (%s) " %(interaction) )

		Base_stack={}

		if(seq_num_1<seq_num_2):
			Base_stack["res1"]=seq_num_1
			Base_stack["res2"]=seq_num_2
			Base_stack["face1"]=face1
			Base_stack["face2"]=face2
		else:
			Base_stack["res1"]=seq_num_2
			Base_stack["res2"]=seq_num_1
			Base_stack["face1"]=face2
			Base_stack["face2"]=face1

		Base_stack["base_ID1"]=SEQUENCE[Base_stack["res1"]-1]	
		Base_stack["base_ID2"]=SEQUENCE[Base_stack["res2"]-1]

		Base_stack["string"]=interaction

		assert_valid_Base_stack(Base_stack)

		Base_stack_list.append(Base_stack)

		if( (double_count) and (seq_num_1 in loop_res_list) and (seq_num_2 in loop_res_list) ):
			Base_stack_list.append( flip_base_stack(Base_stack) )


	while(True):

		line_num+=1

		if(line_num>=len(MC_lines)): break

		line=MC_lines[line_num]

		line_list=line.split()

		if(line.count("--")>0): error_exit_with_message("line.count(\"-\")>0 for line=%s" %(line) )

		seq_num_1=int(line_list[0].split("-")[0][1:])
		seq_num_2=int(line_list[0].split("-")[1][1:])

		if( (seq_num_1 not in loop_res_list) and (seq_num_2 not in loop_res_list) ): continue # not an interaction involving the loop!

		if(seq_num_1==seq_num_2): error_exit_with_message("seq_num_1=(%s)==(%s)=seq_num_2!", (seq_num_1,seq_num_2) )


		Base_pair={}
	
		interaction=line_list[3:]

		num_valid_edge=0

		while(True):

			if(interaction[0].count("/")==0): break 

			potential_edge_1=interaction[0].split("/")[0]
			potential_edge_2=interaction[0].split("/")[1]

			interaction=interaction[1:]

			potential_edge_1=convert_MC_to_LW_base_edge(potential_edge_1, include_bifurcation_BP)
			potential_edge_2=convert_MC_to_LW_base_edge(potential_edge_2, include_bifurcation_BP)

			if(potential_edge_1=="N"): continue
			if(potential_edge_2=="N"): continue

			num_valid_edge+=1

			edge_1=potential_edge_1
			edge_2=potential_edge_2

		if(num_valid_edge==0): continue

		if(num_valid_edge>1): error_exit_with_message("num_valid_edge>1 for line=%s" %(line))

		if(interaction[0]=="adjacent_5p"): interaction=interaction[1:] #Hacky!!
		if(interaction[0]=="inward"):      interaction=interaction[1:] #Hacky!!

		if(len(interaction)<3): error_exit_with_message("len(interaction)=(%d)<3 | interaction=%s | line=%s" %(len(interaction), interaction, line) )

		if(seq_num_1<seq_num_2):
			Base_pair["res1"]=seq_num_1	
			Base_pair["res2"]=seq_num_2	
			Base_pair["edge1"]=edge_1
			Base_pair["edge2"]=edge_2			
		else:
			Base_pair["res1"]=seq_num_2	
			Base_pair["res2"]=seq_num_1	
			Base_pair["edge1"]=edge_2
			Base_pair["edge2"]=edge_1

		Base_pair["base_ID1"]=SEQUENCE[Base_pair["res1"]-1]	
		Base_pair["base_ID2"]=SEQUENCE[Base_pair["res2"]-1]

		if(interaction[2]=="cis"):
			Base_pair["orientation"]="c"
		elif(interaction[2]=="trans"):
			Base_pair["orientation"]="t"
		else:
			error_exit_with_message("Invalid orientation=%s | interaction=%s | line=%s" %(interaction[2], interaction, line) )

		Base_pair["string"]= list_to_string(line_list[2:], first_separator=False)

		assert_valid_Base_pair(Base_pair)

		Base_pair_list.append(Base_pair)

		if( (double_count) and (seq_num_1 in loop_res_list) and (seq_num_2 in loop_res_list) ):
			Base_pair_list.append( flip_base_pair(Base_pair) )



	output_BS_and_BS_list(Base_pair_list, Base_stack_list, title=foldername, output_filename="BP_and_BS_data.txt", To_terminal=Verbose)

	os.chdir( "../")

	return (Base_pair_list, Base_stack_list)

