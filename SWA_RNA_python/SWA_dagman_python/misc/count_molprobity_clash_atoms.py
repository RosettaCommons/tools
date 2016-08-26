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
import string
from os import popen 
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################



#count_molprobity_clash_atoms.py -clash_file CLASH_LIST_STD_BOND_ModeRNA_INPUT_APRIL_22_TLR_C7_2_RLOOM_MODELFHv23.txt > ANALYZE_CLASH_LIST_STD_BOND_ModeRNA_INPUT_APRIL_22_TLR_C7_2_RLOOM_MODELFHv23.txt

#count_molprobity_clash_atoms.py -clash_file CLASH_LIST_STD_BOND_ModeRNA_INPUT_April_22_TLR_ModeRNA_C7_2_modelFHv23.txt > ANALYZE_CLASH_LIST_STD_BOND_ModeRNA_INPUT_April_22_TLR_ModeRNA_C7_2_modelFHv23.txt



#count_molprobity_clash_atoms.py -clash_file CLASH_LIST_ModeRNA_INPUT_APRIL_22_TLR_C7_2_RLOOM_MODELFH.txt > ANALYZE_CLASH_LIST_ModeRNA_INPUT_APRIL_22_TLR_C7_2_RLOOM_MODELFH.txt

#count_molprobity_clash_atoms.py -clash_file Clash_list_ModeRNA_INPUT_April_22_TLR_ModeRNA_C7_2_modelFH.txt > ANALYZED_Clash_list_ModeRNA_INPUT_April_22_TLR_ModeRNA_C7_2_modelFH.txt


#count_molprobity_clash_atoms.py -clash_file MOPROBITY_CLASH_LIST_March_24_3_NT_MUTATTED_RLOOM_MODEL_1EBQ_A_22_26.txt > ANALYZED_MOPROBITY_CLASH_LIST_March_24_3_NT_MUTATTED_RLOOM_MODEL_1EBQ_A_22_26.txt

#count_molprobity_clash_atoms.py -clash_file MOPROBITY_SWA_S_0_000004.txt  > ANALYZED_MOPROBITY_SWA_S_0_000004.txt

#count_molprobity_clash_atoms.py -clash_file MOLPROBITY_CLASH_LIST_March_24_ModeRNA_C7_2_using_2R8S_template.txt > ANALYZED_MOLPROBITY_CLASH_LIST_March_24_ModeRNA_C7_2_using_2R8S_template.txt

#count_molprobity_clash_atoms.py -clash_file TRAIL_2_MOLPROBITY_CLASH_LIST_March_24_MUTATED_RLOOM_MODEL_2GYB_A_649_654.txt  > ANALYZED_TRAIL_2_MOLPROBITY_CLASH_LIST_March_24_MUTATED_RLOOM_MODEL_2GYB_A_649_654.txt

clash_file= parse_options( argv, "clash_file", "" )

if(clash_file==""): error_exit_with_message("clash_file==\"\"!")

if(exists(clash_file)==False): error_exit_with_message("clash_file (%s) doesn't exist!" %(torsion_filename))

include_hydrogen=True

loop_seq_num=[11,12,13]

seq_num_col_list=[3,8]
atom_col_list=[5,10]


num_clash_atoms_list=[]
loop_loop_num_clash_atoms_list=[]


for n in range(0,len(seq_num_col_list)):

	infile=open(clash_file,"r")

	num_clash_atoms=0
	loop_loop_num_clash_atoms=0

	for line in infile:

		if(line=='\n'): continue
		if(line[0]=='#'): continue

		line_split=line.split()
		seq_num=int(line_split[seq_num_col_list[n]])

		if(seq_num not in loop_seq_num): continue



		atom_name_one=line_split[atom_col_list[0]]
		atom_name_two=line_split[atom_col_list[1]]

		if(include_hydrogen==False):
			if(atom_name_one.count("H")!=0): continue

			if(atom_name_two.count("H")!=0): continue


		#Ignore loop-loop interactions.
		seq_num_one=int(line_split[seq_num_col_list[0]])
		seq_num_two=int(line_split[seq_num_col_list[1]])
		if((seq_num_one in loop_seq_num) and (seq_num_two in loop_seq_num) ): 
			loop_loop_num_clash_atoms+=1
		else:
			num_clash_atoms+=1
		
		print line[:-1]


	num_clash_atoms_list.append(num_clash_atoms)
	loop_loop_num_clash_atoms_list.append(loop_loop_num_clash_atoms)

	infile.close()



for n in range(0,len(seq_num_col_list)):
	print "num_atom_clashes method %d=%d" %(n+1, num_clash_atoms_list[n])
	print "loop_loop_num_atom_clashes method %d=%d" %(n+1, loop_loop_num_clash_atoms_list[n])
#if(num_clash_atoms_list[0]!=num_clash_atoms_list[1]): error_exit_with_message("num_clash_atoms_list[0]!=num_clash_atoms_list[1]")



