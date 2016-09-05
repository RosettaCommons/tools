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
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################

def extract_FR3D_data_func(rosetta_pdb_name, sample_segment, double_count=True, allow_near_match=False, Verbose=True):

	matlab="/Applications/MATLAB_R2009b.app/bin/matlab"

	if( exists( rosetta_pdb_name )==False):  error_exit_with_message("rosetta_pdb_name (%s) doesn't exist!" %(rosetta_pdb_name) )

	rosetta_pdb_name=os.path.abspath(rosetta_pdb_name)

	###############################Consistency check#############################################
	fasta_file = "fasta"

	submit_subprocess("pdb2fasta.py %s > %s" %(rosetta_pdb_name, fasta_file) )

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


	###########################consistency check#######################################


	foldername="extract_FR3D_data_%s_%4s_%4s" %( basename(rosetta_pdb_name).replace(".pdb",""), str(start_res).zfill(4), str(end_res).zfill(4))

	if( exists(foldername) ):   
		print "WARNING foldername (%s) already exist! ...removing!..." %(foldername) 
		submit_subprocess("rm -r %s " %(foldername) )

	submit_subprocess( 'mkdir %s' %(foldername)  )

	os.chdir( foldername)


	####################Convert Rosetta format PDB to standard format PDB##########################

	standard_pdb_name="STANDARD_%s" %(basename(rosetta_pdb_name))

	submit_subprocess("Rosetta_to_standard_PDB.py -s %s -output_pdb %s" %(rosetta_pdb_name, standard_pdb_name) )


	###############################################################################################

	sample_res=''
	loop_res_list=[]

	for seq_num in range(start_res, end_res+1):
		sample_res+='%d ' %(seq_num)

		loop_res_list.append(seq_num)

	###############################################################################################
	FR3D_command="%s -nodesktop -nosplash  -r \"extract_BP_and_BS_data('%s')\"" %(matlab, standard_pdb_name)
	print "FR3D_command=%s " %(FR3D_command)
	submit_subprocess( FR3D_command ) #Choose allow sending of multiple pdb?


	###############################################################################################


	FR3D_BP_and_BS_table="BP_and_BS_" + standard_pdb_name[0:-4] + ".txt"

	if(exists(FR3D_BP_and_BS_table)==False): error_exit_with_message("FR3D_BP_and_BS_table (%s) doesn't exist!" %(FR3D_BP_and_BS_table) )

	FR3D_lines=open(FR3D_BP_and_BS_table).readlines()


	FR3D_INFO_line=FR3D_lines[0]

	print "FR3D_INFO_line=", FR3D_INFO_line

	FR3D_SEQ_line=FR3D_lines[1].split()

	FR3D_SEQ=list_to_string( map( lambda x : x[0] , FR3D_SEQ_line ) , separator="")

	print "FR3D_SEQ=" , FR3D_SEQ
	print "SEQUENCE=" , SEQUENCE

	if(len(SEQUENCE)!=len(FR3D_SEQ_line)): error_exit_with_message("len(SEQUENCE)=(%s)!=(%s)=len(FR3D_SEQ_line))" %(len(SEQUENCE),len(FR3D_SEQ_line)) )

	if(len(SEQUENCE)!=len(FR3D_SEQ)): error_exit_with_message("len(SEQUENCE)=(%s)!=(%s)=len(FR3D_SEQ))" %(len(SEQUENCE),len(FR3D_SEQ)) )



	for n in range(len(FR3D_SEQ_line)):
		if(SEQUENCE[n]!=FR3D_SEQ[n]): error_exit_with_message("SEQEUNCE[n]=(%s)!=(%s)=FR3D_SEQ_line[n]" %(SEQEUNCE[n],FR3D_SEQ_line[n]) )

		seq_num=n+1

		if(seq_num!=int(FR3D_SEQ_line[n][1:])): error_exit_with_message("seq_num=(%s)!=(%s)=int(FR3D_SEQ_line[n][1:]" %(seq_num,int(FR3D_SEQ_line[n][1:]) ) )


	Base_pair_list=[]
	Base_stack_list=[]

	for n in range(1, len(FR3D_lines)):

		print FR3D_lines[n],

	for n in range(2, len(FR3D_lines)):

		line=FR3D_lines[n]
		line_list=line.split()

		seq_num_1=n-1

		if(line_list[0][0]!=SEQUENCE[seq_num_1-1]): error_exit_with_message("line_list[0][0]=(%s)!=(%s)=SEQUENCE[seq_num_1-1]" %(line_list[0][0], SEQUENCE[seq_num_1-1]) )
	
		if(line_list[seq_num_1][0]!=SEQUENCE[seq_num_1-1]): error_exit_with_message("line_list[seq_num_1-1][0]=(%s)!=(%s)=SEQUENCE[seq_num_1-1]" %(line_list[seq_num_1-1][0], SEQUENCE[seq_num_1-1]) )

		if(len(line_list)!=len(SEQUENCE)+1): error_exit_with_message("len(line_list)=(%s)!=(%s)=len(SEQUENCE)+1" %(len(line_list),len(SEQUENCE)+1) ) 

		for ii in range( seq_num_1+1, len(SEQUENCE)+1 ):

			seq_num_2=ii

			if( (seq_num_1 not in loop_res_list) and (seq_num_2 not in loop_res_list) ): continue # not an interaction involving the loop!

			if(seq_num_1==seq_num_2): error_exit_with_message("seq_num_1=(%s)==(%s)=seq_num_2!", (seq_num_1,seq_num_2) )


			interaction=line_list[seq_num_2]

			if(interaction=="----"): continue #no interaction!

			if(interaction[0]=='n'):
				if(allow_near_match):
					interaction=interaction[1:]
				else:
		 			 continue #ignore near-match for now

			if(interaction[0:3]=="bif"): continue #ignore bifurcation!

			if(interaction[0]=='s'): #Base_stack!
				Base_stack={}
				if(seq_num_1 in loop_res_list):
					Base_stack["res1"]=seq_num_1		
					Base_stack["res2"]=seq_num_2	
					Base_stack["face1"]=interaction[1]
					Base_stack["face2"]=interaction[2]
				else:
					Base_stack["res1"]=seq_num_2		
					Base_stack["res2"]=seq_num_1	
					Base_stack["face1"]=interaction[2]
					Base_stack["face2"]=interaction[1]

				Base_stack["base_ID1"]=SEQUENCE[Base_stack["res1"]-1]	
				Base_stack["base_ID2"]=SEQUENCE[Base_stack["res2"]-1]

				Base_stack["string"]=line_list[seq_num_2]

				assert_valid_Base_stack(Base_stack)

				Base_stack_list.append(Base_stack)

				if( (double_count) and (seq_num_1 in loop_res_list) and (seq_num_2 in loop_res_list) ):
					Base_stack_list.append( flip_base_stack(Base_stack) )

			elif(interaction[0]=='t' or interaction[0]=='c'): #Base_pair!
				Base_pair={}

				if(seq_num_1 > seq_num_2):
					Base_pair["res1"]=seq_num_1	
					Base_pair["res2"]=seq_num_2	
					Base_pair["edge1"]=interaction[1].upper()
					Base_pair["edge2"]=interaction[2].upper()
				else:
					Base_pair["res1"]=seq_num_2	
					Base_pair["res2"]=seq_num_1	
					Base_pair["edge1"]=interaction[2].upper()
					Base_pair["edge2"]=interaction[1].upper()

				Base_pair["base_ID1"]=SEQUENCE[Base_pair["res1"]-1]	
				Base_pair["base_ID2"]=SEQUENCE[Base_pair["res2"]-1]

				Base_pair["orientation"]=interaction[0]
				Base_pair["string"]=line_list[seq_num_2]

				assert_valid_Base_pair(Base_pair)

				Base_pair_list.append(Base_pair)

				if( (double_count) and (seq_num_1 in loop_res_list) and (seq_num_2 in loop_res_list) ):
					Base_pair_list.append( flip_base_pair(Base_pair) )

			else:
				error_exit_with_message("Invalid interaction (%s)!, seq_num_1=%d, seq_num_2=%d" %(interaction, seq_num_1, seq_num_2) )


	output_BS_and_BS_list(Base_pair_list, Base_stack_list, title=foldername, output_filename="BP_and_BS_data.txt", To_terminal=Verbose)

	return (Base_pair_list, Base_stack_list)

