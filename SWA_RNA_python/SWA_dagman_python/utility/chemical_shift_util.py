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
##################################################################

from SWA_util import *
import PDB_operations 

##################################################################
def extract_experimental_chemical_shift(input_BMRB, proton_name_list, SEQUENCE):

	ALL_exp_shift_data_list=[]

	BMRB_LINE_LIST=open( input_BMRB ).readlines() 
	
	for line in BMRB_LINE_LIST:

		line_list=line.split()

		if(len(line_list)!=9): error_exit_with_message("len(line_list)=(%s)!=9 for line=(%s)" %(len(line_list), line) )

		shift_data={}

		shift_data["rsd_name"]  =line_list[3]		

		if(line_list[1]!=line_list[2]): error_exit_with_message("line_list[1]!=line_list[2] for line=%s" %(line))

		BMRB_seq_num=int(line_list[1])

		#if(BMRB_to_PDB_seq_num_map.has_key(BMRB_seq_num)==False): continue

		shift_data["seq_num"]   =BMRB_seq_num

		shift_data["atom_name"] =line_list[4]
		shift_data["calc_shift"]="BLAH"
		shift_data["exp_shift"] =float(line_list[6])

		if(shift_data["rsd_name"]!=SEQUENCE[shift_data["seq_num"]-1]): 
			error_exit_with_message("shift_data[\"rsd_name\"]=(%s)!=(%s)=SEQUENCE[%d-1]" %(shift_data["rsd_name"],  SEQUENCE[shift_data["seq_num"]-1], shift_data["seq_num"] ) )

		if(shift_data["atom_name"].count('H')>0):#If hydrogen atom
			PDB_operations.assert_is_valid_standard_hatom_name(shift_data["atom_name"], matching_whitespace=False)

		if(proton_name_list.count(shift_data["atom_name"])==0): continue

		ALL_exp_shift_data_list.append(shift_data)

	return ALL_exp_shift_data_list

