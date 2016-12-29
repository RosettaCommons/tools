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
########################################################################

from SWA_dagman_python.utility.SWA_util import *
########################################################################


########################################################################

def Is_nested_inside(chunk_info_one, chunk_info_two): #check whether chunk_info_one is nested within chunk_info_two

	if(chunk_info_one["res_list"][0]>chunk_info_two["res_list"][0] and chunk_info_one["res_list"][-1]<chunk_info_two["res_list"][-1]):
		return True
	else:
		return False 


########################################################################
def check_valid_bracket_structure(dot_bracket_structure, total_res):

	first_strand_list=[]
	second_strand_list=[]

	num_open_brackets=0
	num_close_brackets=0

	if(len(dot_bracket_structure)!=total_res): error_exit_with_message("len(dot_bracket_structure)!=total_res")

	for n in range(total_res):
		seq_num=n+1
		char=dot_bracket_structure[n]

		if(char=="(" and (seq_num==1 or dot_bracket_structure[n-1]!="(")):
			first_strand_list.append([])

		if(char=="("):
			first_strand_list[-1].append(seq_num)
			num_open_brackets+=1

	for n in range(total_res-1,-1,-1):
		seq_num=n+1		
		char=dot_bracket_structure[n]

		if(char==")" and (seq_num==total_res or dot_bracket_structure[n+1]!=")")):
			second_strand_list.append([])

		if(char==")"):
			second_strand_list[-1].append(seq_num)
			num_close_brackets+=1




	if(len(first_strand_list)!=len(second_strand_list)):
		error_exit_with_message("len(first_strand_list)=(%s)!=(%s)=len(second_strand_list)" %(first_strand_list, second_strand_list) )

	if(num_open_brackets!=num_close_brackets): error_exit_with_message("num_open_brackets(%d)!=num_close_brackets(%d)" %(num_open_brackets, num_close_brackets) )

########################################################################


def parse_dot_bracket_structure(dot_bracket_structure, total_res):

	dot_bracket_structure=dot_bracket_structure.replace("[","(")
	dot_bracket_structure=dot_bracket_structure.replace("]",")")

	if(len(dot_bracket_structure)!=total_res):
		error_exit_with_message("len(dot_bracket_structure) !=total_res (%d), len(dot_bracket_structure)=%s" %(total_res, len(dot_bracket_structure)) )

	check_valid_bracket_structure(dot_bracket_structure, total_res)


	bracket_pair_list=[]

	accounted_pos=[]

	remaining_dot_bracket_structure=dot_bracket_structure



	while(("(" in remaining_dot_bracket_structure) or (")" in remaining_dot_bracket_structure) ):

		if(len(dot_bracket_structure)!=total_res): error_exit_with_message("len(dot_bracket_structure)!=total_res")

		curr_open_bracket=[]
		curr_close_bracket=[]

		for n in range(total_res):
			seq_num=n+1
			char=dot_bracket_structure[n]
			if(n in accounted_pos): continue
			
			if(char=="(" and (seq_num==1 or dot_bracket_structure[n-1]!="(")): #start of a new open_bracket!
				curr_open_bracket=[]

			if(char=="("):
				curr_open_bracket.append(seq_num)

			if(char==")"):
				curr_close_bracket.append(seq_num)

			if(char==")" and (seq_num==total_res or dot_bracket_structure[n+1]!=")")): #end of the close_bracket..this should match the last open_bracket
				first_open=curr_open_bracket[0]-1
				last_close=curr_close_bracket[-1]-1

				if(len(curr_open_bracket)!=len(curr_close_bracket)): 
					error_exit_with_message("len(curr_open_bracket)=(%d)!=(%d)=len(curr_close_bracket), location: ~%d and ~%d" %(len(curr_open_bracket), len(curr_close_bracket), first_open, last_close) )
					
				bracket_pair_list.append([curr_open_bracket,curr_close_bracket]	)
				accounted_pos.extend(range( first_open, last_close+1 ) )
				
				break

		remaining_dot_bracket_structure=""
		for n in range(len(dot_bracket_structure)):
			if(n not in accounted_pos):
				remaining_dot_bracket_structure=remaining_dot_bracket_structure+dot_bracket_structure[n]

	chunk_res_list=[]


	for n in range(len(bracket_pair_list)):
		chunk_res=[]
		chunk_res.extend(bracket_pair_list[n][0])
		chunk_res.extend(bracket_pair_list[n][1])
		chunk_res.sort()
		chunk_res_list.append(chunk_res)

	chunk_res_list.sort()

	print "bracket_pair_list=", bracket_pair_list
	#print "chunk_res_list=", chunk_res_list
	

	return chunk_res_list

########################################################################

def get_chunk_res_list_from_cutpoint_open(cutpoint_open_list, helix_length_list, total_res):

	if(len(helix_length_list)!=0):
		for n in range(len(helix_length_list)): 
			if(helix_length_list[n]<0): error_exit_with_message("ERROR! helix_length_list[%d]<0!" %(n) )
		if((len(cutpoint_open_list)+1)!=len(helix_length_list)):
			error_exit_with_message("len(cutpoint_open_list)+1)=%s!=%s=len(helix_length_list)" %(len(cutpoint_open_list)+1, len(helix_length_list) ) )
	else:
		for n in range(len(cutpoint_open_list)+1): 
			helix_length_list.append(2)

	chunk_res_list=[]

	for n in range(len(helix_length_list)):
		helix_length=helix_length_list[n]

		chunk_res=[]
		if(n==0):
			strand_1_seq_num=range(1,helix_length+1)
			strand_2_seq_num=range(total_res-helix_length+1,total_res+1)	
		else:
			cutpoint_open=cutpoint_open_list[n-1]
			strand_1_seq_num=range(cutpoint_open-helix_length+1,cutpoint_open+1)
			strand_2_seq_num=range(cutpoint_open+1,cutpoint_open+helix_length+1)

		chunk_res.extend(strand_1_seq_num)
		chunk_res.extend(strand_2_seq_num)	
		chunk_res_list.append(chunk_res)

	return chunk_res_list


########################################################################
def get_sample_regions(chunk_res_list, virtual_phosphate_list, total_res):

	prev_in_chunk=[]
	curr_in_chunk=[]
	sample_res_list=[]
	sample_bb_list=[] #the backbone 5' of the residue sugar.

	for seq_num in range(1, total_res+1):

		for n in range(len(chunk_res_list)):
			if(seq_num in chunk_res_list[n]):
				curr_in_chunk.append(n)

		if(len(curr_in_chunk)==0): 
			sample_res_list.append(seq_num)
			sample_bb_list.append(seq_num)
	
		if	(len(curr_in_chunk)!=0):
			all_chunks_are_discontinuous=True
			for chunk_ID in curr_in_chunk:
				if(chunk_ID in prev_in_chunk):
					all_chunks_are_discontinuous=False

			if(all_chunks_are_discontinuous): sample_bb_list.append(seq_num)

		prev_in_chunk=curr_in_chunk
		curr_in_chunk=[]	

	if(len(sample_res_list)==0): error_exit_with_message("len(sample_res_list)==0")
	if(len(sample_bb_list)==0): error_exit_with_message("len(sample_bb_list)==0")
	
	#Don't sample these bb torsions since they will be virtualized anyways!
	sample_bb_list=list( Set(sample_bb_list)-Set(virtual_phosphate_list) )
	sample_bb_list.sort()

	return (sample_res_list, sample_bb_list)

########################################################################

def get_allow_bulge_res_list(allow_bulge_res_list, sample_res_list, cutpoint_open_list, total_res):

	#Not sure why not desirable to add virtual_res to 5' and 3' of cutpoint_open and final_res. Might want to allow this in the FUTURE!
	if(len(allow_bulge_res_list)==0): #if allow didn't user specify
		allow_bulge_res_list=list(Set(sample_res_list)-Set(cutpoint_open_list)-Set(map( lambda x : x+1 , cutpoint_open_list ))-Set([total_res]))
		allow_bulge_res_list.sort()

	if(len(allow_bulge_res_list)==0): error_exit_with_message("len(allow_bulge_res_list)==0")

	return allow_bulge_res_list

########################################################################
def get_res_list(segment_string):

	if(len(segment_string.split("-"))!=2):error_exit_with_message("len(segment_string.split(\"-\"))!=2, segment_string=%s" %(segment_string) )

	start_res=int(segment_string.split("-")[0])

	end_res=int(segment_string.split("-")[1])

	if(start_res>end_res): error_exit_with_message("start_res(%d)>end_res(%d)" %(start_res, end_res))

	res_list=[]

	for seq_num in range(start_res, end_res+1):
		res_list.append(seq_num)	

	return res_list


def get_chunk_files(chunk_info, input_filename, chunk_ID, sequence, sample_res_list, cutpoint_open_list):

		if(exists("SETUP_LOG/")==False): submit_subprocess("mkdir SETUP_LOG/")

		VDW_screener_info=[]

		##############################################################################################################################
		if(input_filename!="IDEAL_HELIX"):
			if(exists(input_filename)==False): error_exit_with_message("input_filename (%s) doesn't exist!" %(input_filename) )
			chunk_info["filename"]=input_filename
		else:
			#create a idealize helix!
			if(len(chunk_info["segment_list"])!=2): 
				print "chunk_info[\"segment_list\"]=", chunk_info["segment_list"]
				error_exit_with_message("len(chunk_info[\"segment_list\"])!=2")
			
			strand_one=get_res_list(chunk_info["segment_list"][0])
			strand_two=get_res_list(chunk_info["segment_list"][1])

			if(len(strand_one)!=len(strand_two)): error_exit_with_message("len(strand_one)=(%d)!=(%d)=len(strand_two)" %(len(strand_one), len(strand_two)) )

			VDW_screener_type="NONE"
			num_valid_VDW_rep_type=0

			if( ((strand_one[0]-1) not in sample_res_list) and ((strand_two[-1]+1) not in sample_res_list) ):
				VDW_screener_type="LOWER"	
				num_valid_VDW_rep_type+=1
						
			if( (strand_one[-1] in cutpoint_open_list) and ( (strand_one[-1]+1)==strand_two[0] ) ):
				VDW_screener_type="UPPER"	
				num_valid_VDW_rep_type+=1

			if(num_valid_VDW_rep_type>1): error_exit_with_message("num_valid_VDW_rep_type (%d) > 1 " %(num_valid_VDW_rep_type) )

			##############################################################
			verbose=True
			helix_filename="helix_stub_%s.pdb" %(chunk_ID)

			create_helix_command="create_idealize_helix_general.py "
			create_helix_command+="-sequence %s -verbose %s -VDW_screener_type %s " %(sequence, verbose, VDW_screener_type)
			create_helix_command+="-strand_1_seq_num %s -strand_2_seq_num %s -helix_filename %s" %(list_to_string(strand_one), list_to_string(strand_two), helix_filename) 
			create_helix_command+="> SETUP_LOG/LOG_create_helix_%d.txt " %(chunk_ID)

			if(exists('VDW_rep_screen_info.txt')==True): submit_subprocess("rm -r VDW_rep_screen_info.txt")

			submit_subprocess(create_helix_command) 

			chunk_info["filename"]=helix_filename

			if(VDW_screener_type!="NONE"):
				VDW_rep_screen_info_lines=open( "VDW_rep_screen_info.txt" ).readlines()
				if(len(VDW_rep_screen_info_lines)!=1): error_exit_with_message("len(VDW_rep_screen_info_lines)!=1")
				VDW_screener_info=VDW_rep_screen_info_lines[0].split()[1:]
			##############################################################

		if(chunk_info["filename"][-4:]==".pdb"): ##--> convert pdb to silent_file!!
			silent_file=chunk_info["filename"][0:-4] + ".out"
			if(exists(silent_file)): error_exit_with_message("silent_file (%s) already exist!" %(silent_file) )

			tag_name=basename(silent_file[0:-4])

			pdb_to_silent_file_command="pdb_to_silent_file.py -pdb %s -tag_name %s -output_silent_file %s " %(chunk_info["filename"], tag_name, silent_file)
			pdb_to_silent_file_command+="> SETUP_LOG/LOG_pdb_to_silent_%d.out" %(chunk_ID)
			submit_subprocess(pdb_to_silent_file_command)

			chunk_info["filename"]=silent_file

		#check that file is a silent_file!
		if(chunk_info["filename"][-4:]!=".out"):
			error_exit_with_message("chunk_info[\"filename\"] should be a silent_file, but doesn't have .out extension!")

		###############################################################################################################################
		return ( chunk_info["filename"], VDW_screener_info )
		#return ( os.path.abspath(chunk_info["filename"]), VDW_screener_info )


########################################################################

def get_chunk_info_list(chunk_res_list, chunk_file_list, sequence, sample_res_list, cutpoint_open_list):
	######################Initial the chunk_info_list##################################
	#OK for each chunk, the following information is needed
	#1. chunk_info["res_list"]
	#2. chunk_info["segment_list"]
	#3. chunk_info["filename"] #either a pdb or a silent_file.
	#4. chunk_info["VDW_screener_info"] #pdb for VDW_rep_screening if chunk_info["filename"] is a idealized helix at cutpoint or outer!

	chunk_info_list=[]

	if(chunk_file_list==[""]): 
		chunk_file_list=	map( lambda x : "IDEAL_HELIX" , chunk_res_list )

	if(len(chunk_file_list)!=len(chunk_res_list)):
		error_exit_with_message("len(chunk_file_list)=(%s)!=(%s)=len(chunk_res_list)" %(len(chunk_file_list), len(chunk_res_list)) )

	for n in range(len(chunk_res_list)):

		chunk_info={}	

		chunk_res_list[n].sort()

		chunk_info["res_list"]=copy.deepcopy(chunk_res_list[n])

		chunk_info["segment_list"]=get_segment_string_list(chunk_info["res_list"], cutpoint_open_list)		

		(chunk_info["filename"],chunk_info["VDW_screener"])=get_chunk_files(chunk_info, chunk_file_list[n], n, sequence, sample_res_list, cutpoint_open_list)   
		
		chunk_info_list.append(chunk_info)

	print "------------------------------chunk_info_list------------------------------"
	for n in range(len(chunk_info_list)):
		print " chunk_info #%d: " %(n), chunk_info_list[n] 
	print "---------------------------------------------------------------------------"

	return chunk_info_list

########################################################################


def get_native_align_res_list(long_loop_mode, sample_res_list, sample_bb_list, total_res, user_native_alignment_res):

	all_res=range(1,total_res+1)

	native_alignment_res=[]

	if(len(user_native_alignment_res)>0):
		native_alignment_res=user_native_alignment_res
		print "USER passed in native_alignment_res %s" %(list_to_string(native_alignment_res))

	elif(long_loop_mode): #prefect alignment to base_pdb in long_loop_mode..
		native_alignment_res=[1, total_res]
	else:
		for seq_num in all_res: #basically the sample_res and its intermediate neighboring res.
			if( (seq_num in sample_res_list)     and (seq_num not in native_alignment_res) ): native_alignment_res.append(seq_num)
			if( ((seq_num+1) in sample_res_list) and (seq_num not in native_alignment_res) ): native_alignment_res.append(seq_num)
			if( ((seq_num-1) in sample_res_list) and (seq_num not in native_alignment_res) ): native_alignment_res.append(seq_num)

			if( (seq_num in sample_bb_list)     and (seq_num not in native_alignment_res) ): native_alignment_res.append(seq_num)
			if( ((seq_num+1) in sample_bb_list) and (seq_num not in native_alignment_res) ): native_alignment_res.append(seq_num)

	native_alignment_res.sort()

	if(len(native_alignment_res)==0): error_exit_with_message("len(native_alignment_res)==0")

	return native_alignment_res

########################################################################
def get_jump_points_within_chunk_segments(chunk_info_list, fixed_res_list):

	jump_point_list=[]
	additional_cutpoints_close=[]

	for chunk_ID in range(len(chunk_info_list)): #Jump within chunk (rare!)

		chunk_info=chunk_info_list[chunk_ID]

		for n in range(len(chunk_info["segment_list"])):

			segment=		map( lambda x : int(x) , chunk_info["segment_list"][n].split("-") )

			add_jump_point=False

			#If the first and last res of the segment is fixed and there are minimize_res in between. Then add jump-point to keep the ends fixed!
			if( (segment[0] in fixed_res_list) and (segment[1] in fixed_res_list) ): 
				for seq_num in range(segment[0], segment[1]+1):
					if(seq_num not in fixed_res_list): add_jump_point=True

			if(add_jump_point==False): continue

			for seq_num in range(segment[0], segment[1]+1):
				if(seq_num not in fixed_res_list):
					if( (seq_num)   not in additional_cutpoints_close): additional_cutpoints_close.append(seq_num)
					if( (seq_num-1) not in additional_cutpoints_close): additional_cutpoints_close.append(seq_num-1)

			print "Add additional jump_point within chunk segment:" , chunk_info, " | segment= ", chunk_info["segment_list"][n]
		
			five_prime_jump =segment[0]
			three_prime_jump=segment[1]

			jump_point_list.append([five_prime_jump, three_prime_jump])

	jump_point_list.sort()		
	additional_cutpoints_close.sort()
	print "jump_points_within_chunk_segments=", jump_point_list, " | additional_cutpoints_close=", additional_cutpoints_close

	return (jump_point_list, additional_cutpoints_close)

########################################################################

def get_jump_points_between_chunk_segments(chunk_info_list, sample_res_list, sample_bb_list, fixed_res_list):

	jump_point_list=[]

	for chunk_ID in range(len(chunk_info_list)): #Jump between chunk!
	
		chunk_info=chunk_info_list[chunk_ID]

		for n in range(len(chunk_info["segment_list"])-1):
			first_segment=		map( lambda x : int(x) , chunk_info["segment_list"][n].split("-") )
			second_segment=	map( lambda x : int(x) , chunk_info["segment_list"][n+1].split("-") )
	
			five_prime_jump=first_segment[0]
			three_prime_jump=second_segment[1] 

			while(True):
				if( five_prime_jump < first_segment[0] ): error_exit_with_message("five_prime_jump < first_segment[0]!")
				if( five_prime_jump > first_segment[1] ): error_exit_with_message("five_prime_jump > first_segment[1]!")

				valid_jump=True
				if( five_prime_jump in sample_res_list ): valid_jump=False
				if( five_prime_jump not in fixed_res_list): valid_jump=False

				if(valid_jump): break
				
				five_prime_jump+=1

			while(True):
				if( three_prime_jump < second_segment[0] ): error_exit_with_message("three_prime_jump < second_segment[0]!")
				if( three_prime_jump > second_segment[1] ): error_exit_with_message("three_prime_jump > second_segment[1]!")

				valid_jump=True
				if( three_prime_jump in sample_res_list ): valid_jump=False
				if( three_prime_jump not in fixed_res_list): valid_jump=False

				if(valid_jump): break
				
				three_prime_jump-=1

			######Hacky June 02, 2011########################################################################################
			#Basically if there is no sample_res between this chunk and the one inner to it..then don't need to add jump_point.
			continuous_with_other_chunk=False
			other_chunk_ID=-1

			for ii in range(len(chunk_info_list)):
				if(ii==chunk_ID): continue
				other_chunk_info=chunk_info_list[ii]

				if( Is_disjoint_list(range(first_segment[0], first_segment[1]+1), other_chunk_info["res_list"])==False ):
					if( Is_disjoint_list(range(second_segment[0], second_segment[1]+1), other_chunk_info["res_list"])==False ):
						if( first_segment[0] <= min(other_chunk_info["res_list"]) and second_segment[1] >= max(other_chunk_info["res_list"]) ):
							continuous_with_other_chunk=True
							other_chunk_ID=ii

			if(continuous_with_other_chunk):	
				print "strands %d and %d of chunk #%d is continuous with chunk #%d, no jump point added!" %(n, n+1, chunk_ID, other_chunk_ID)
				continue 
			#############################################################################

			jump_point_list.append([five_prime_jump, three_prime_jump])

	jump_point_list.sort()		
	print "jump_points_between_chunk_segments=", jump_point_list

	return jump_point_list

########################################################################

def get_fixed_res_list(long_loop_mode, fixed_res_list, sample_res_list, sample_bb_list, total_res):
#####--> Eventually change this to minimize_res_list and minimize_bb_list and pass into params file.

	minimize_res_list=[]

	all_res=range(1,total_res+1)

	if(len(fixed_res_list)==0): #user did not pass in fixed_res
		
		for seq_num in all_res:
			
			if(seq_num in sample_res_list): minimize_res_list.append(seq_num)
	
			if(long_loop_mode==False):
				if( ((seq_num-1) in sample_res_list) and (seq_num not in minimize_res_list) ): minimize_res_list.append(seq_num)
				if( ((seq_num+1) in sample_res_list) and (seq_num not in minimize_res_list) ): minimize_res_list.append(seq_num)

		fixed_res_list=list(Set(all_res)-Set(minimize_res_list))
		fixed_res_list.sort() 
	else:
		minimize_res_list=list(Set(all_res)-Set(fixed_res_list))
		minimize_res_list.sort() 		

	print_seq_num_list("minimize_res_list=", minimize_res_list)
	print_seq_num_list("fixed_res_list=", fixed_res_list)

	return fixed_res_list

########################################################################

