#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
import os
import copy
#######################


from SWA_util import *
#######################


'''
LP:lone pair pseudo-atoms 

Base_pair={}
Base_pair["res1"]
Base_pair["res2"]
Base_pair["edge1"] ("S" sugar, "W" watson-crick, "H" hoogsteen)
Base_pair["edge2"] ("S" sugar, "W" watson-crick, "H" hoogsteen)
Base_pair["orientation"] //Leontis Westhof base-pair orientation ("c": cis, "t": trans) 
Base_pair["string"] //FR3D string

Base_stack={}
Base_stack["res1"]
Base_stack["res2"]
Base_stack["face1"] ("5", "3")
Base_stack["face2"] ("5", "3")
Base_stack["string"] //FR3D string

Convert between MC_annotate and FR3D Base-stack orientations: upward==s35, downward==s53, inward==s33, outward==s55

Type s35 for stacking in which the first base uses its 3 face, and the second base uses its 5 face. Similarly, type s53, s33, or s55. Type "stack" to allow all stacking interactions. 
'''
RNA_BP_AND_BS_UTIL_ALLOW_BIFUCATION=True
########################################################################

def 	Is_valid_base_ID(base_ID):

	if(base_ID=="A"): return True
	if(base_ID=="G"): return True
	if(base_ID=="C"): return True
	if(base_ID=="U"): return True

	return False

def Is_valid_base_edge(base_edge):

	if(base_edge=="W"): return True

	if(base_edge=="H"): return True

	if(base_edge=="S"): return True

	if(RNA_BP_AND_BS_UTIL_ALLOW_BIFUCATION):
		if(base_edge=="Bs"): return True
		if(base_edge=="Bh"): return True

	return False


def assert_valid_base_ID(base_ID):

	if(Is_valid_base_ID(base_ID)==False): error_exit_with_message("Invalid base_ID (%s) " %(base_ID) )

########################################################################

def assert_valid_Base_stack(Base_stack):

	if(Base_stack.has_key("res1")==False): 
		print "Problem Base_stack=", Base_stack
		error_exit_with_message('Base_stack.has_key("res1")==False')

	if(Base_stack.has_key("res2")==False): 
		print "Problem Base_stack=", Base_stack
		error_exit_with_message('Base_stack.has_key("res2")==False')

	if(Base_stack["res1"]==Base_stack["res2"]): error_exit_with_message("Base_stack[\"res1\"]=(%s)==(%s)=Base_stack[\"res2\"]", (Base_stack["res1"], Base_stack["res2"]) )

	if(Base_stack["face1"]!="5" and Base_stack["face1"]!="3"): 
		print "Problem Base_stack=", Base_stack
		error_exit_with_message('Base_stack["face1"]!="5" and Base_stack["face1"]!="3"')
	
	if(Base_stack["face2"]!="5" and Base_stack["face2"]!="3"): 
		print "Problem Base_stack=", Base_stack
		error_exit_with_message('Base_stack["face2"]!="5" and Base_stack["face1"]2="3"')

	if(Is_valid_base_ID(Base_stack["base_ID1"])==False):
		print "Problem Base_pair=", Base_stack
		error_exit_with_message('Invalid Base_stack["Base_ID"]=%s' %(Base_stack["Base_ID"]))

	if(Is_valid_base_ID(Base_stack["base_ID2"])==False):
		print "Problem Base_pair=", BBase_stack
		error_exit_with_message('Invalid Base_stack["Base_ID"]=%s' %(Base_stack["Base_ID"]))

	if(Is_equivalent_base_stack(Base_stack, Base_stack, assert_valid=False )==False): 
		print "Problem Base_stack=", Base_stack
		error_exit_with_message("Base_stack is not equivalent to itself!")

########################################################################


def assert_valid_Base_pair(Base_pair):

	if(Base_pair.has_key("res1")==False): 
		print "Problem Base_pair=", Base_pair 
		error_exit_with_message('Base_pair.has_key("res1")==False')

	if(Base_pair.has_key("res2")==False): 
		print "Problem Base_pair=", Base_pair 
		error_exit_with_message('Base_pair.has_key("res2")==False')

	if(Base_pair["res1"]==Base_pair["res2"]): error_exit_with_message("Base_pair[\"res1\"]=(%s)==(%s)=Base_pair[\"res2\"]", (Base_pair["res1"], Base_pair["res2"]) )

	if(Is_valid_base_edge(Base_pair["edge1"])==False):	
		print "Problem Base_pair=", Base_pair 
		error_exit_with_message('Invalid Base_pair["edge1"]=%s' %(Base_pair["edge1"]))

	if(Is_valid_base_edge(Base_pair["edge2"])==False):	
		print "Problem Base_pair=", Base_pair 
		error_exit_with_message('Invalid Base_pair["edge1"]=%s' %(Base_pair["edge1"]))

	if(Is_valid_base_ID(Base_pair["base_ID1"])==False):
		print "Problem Base_pair=", Base_pair 
		error_exit_with_message('Invalid Base_pair["Base_ID1"]=%s' %(Base_pair["Base_ID1"]))

	if(Is_valid_base_ID(Base_pair["base_ID2"])==False):
		print "Problem Base_pair=", Base_pair 
		error_exit_with_message('Invalid Base_pair["Base_ID2"]=%s' %(Base_pair["Base_ID2"]))

	if(Base_pair["orientation"]!="t" and Base_pair["orientation"]!="c"):
		print "Problem Base_pair=", Base_pair 
		error_exit_with_message('Base_pair["orientation"]!="t" and Base_pair["orientation"]!="c"')
	
	if(Is_equivalent_base_pair(Base_pair, Base_pair, assert_valid=False)==False): 
		print "Problem Base_pair=", Base_pair
		error_exit_with_message("Base_pair is not equivalent to itself!")

########################################################################


def flip_base_pair(start_base_pair):

	flipped_base_pair={}
	flipped_base_pair["res1"]					=start_base_pair["res2"]
	flipped_base_pair["res2"]					=start_base_pair["res1"]
	flipped_base_pair["edge1"]				=start_base_pair["edge2"]
	flipped_base_pair["edge2"]				=start_base_pair["edge1"]
	flipped_base_pair["base_ID1"]			=start_base_pair["base_ID2"]
	flipped_base_pair["base_ID2"]			=start_base_pair["base_ID1"]
	flipped_base_pair["orientation"]	=start_base_pair["orientation"]
	flipped_base_pair["string"]				=start_base_pair["string"]

	assert_valid_Base_pair(flipped_base_pair)

	return flipped_base_pair

########################################################################

def flip_base_stack(start_base_stack):

	flipped_base_stack={}
	flipped_base_stack["res1"]			=start_base_stack["res2"]
	flipped_base_stack["res2"]			=start_base_stack["res1"]
	flipped_base_stack["face1"]			=start_base_stack["face2"]
	flipped_base_stack["face2"]			=start_base_stack["face1"]
	flipped_base_stack["base_ID1"]	=start_base_stack["base_ID2"]
	flipped_base_stack["base_ID2"]	=start_base_stack["base_ID1"]
	flipped_base_stack["string"]		=start_base_stack["string"]

	assert_valid_Base_stack(flipped_base_stack)

	return flipped_base_stack

########################################################################
def output_BS_and_BS_list(Base_pair_list, Base_stack_list, title, output_filename, To_terminal=True):

	Base_pair_list = sort_dictionary_list(Base_pair_list, "res1") 
	Base_stack_list= sort_dictionary_list(Base_stack_list, "res1") 

	BP_and_BS_DATA = open( output_filename, 'w')

	BP_and_BS_DATA.write( "----------------------%s----------------------\n" %(title) )

	for n in range(len(Base_stack_list)):
		#print "Base_stack #%d=" %(n+1), Base_stack_list[n]

		#Warning this doesn't work in python 2.3!
		BP_and_BS_DATA.write( "Base_stack #%2d {" %(n+1) )
		for key in sorted(Base_stack_list[n].iterkeys()):
			BP_and_BS_DATA.write( "'%4s': '%4s'" % (key, Base_stack_list[n][key]) )
		BP_and_BS_DATA.write( "}\n")



	for n in range(len(Base_pair_list)):
		#print "Base_pair #%d=" %(n+1), Base_pair_list[n]

		#Warning this doesn't work in python 2.3!
		BP_and_BS_DATA.write( "Base_pair #%2d {" %(n+1) )
		for key in sorted(Base_pair_list[n].iterkeys()):
			BP_and_BS_DATA.write( "'%4s': '%4s'" % (key, Base_pair_list[n][key]) )
		BP_and_BS_DATA.write( "}\n")

	BP_and_BS_DATA.write( "----------------------" )

	for n in range(len(title)):
		BP_and_BS_DATA.write( "-" )

	BP_and_BS_DATA.write( "----------------------\n" )

	BP_and_BS_DATA.close()

	if(To_terminal):
		bp_and_bs_lines = open( "BP_and_BS_data.txt", 'r').readlines()

		for line in bp_and_bs_lines:

			print line,

########################################################################

def Is_equivalent_base_stack(base_stack_1, base_stack_2, assert_valid=True):

	if(assert_valid):
		assert_valid_Base_stack(base_stack_1)
		assert_valid_Base_stack(base_stack_2)

	return ( Is_equivalent_base_stack_strict(base_stack_1, base_stack_2) ) #Does not treat flip basepair as equivalent!

########################################################################

def Is_equivalent_base_stack_strict(base_stack_1, base_stack_2):

	if(base_stack_1["res1"]!=base_stack_2["res1"]): return False

	if(base_stack_1["res2"]!=base_stack_2["res2"]): return False

	if(base_stack_1["base_ID1"]!=base_stack_2["base_ID1"]): error_exit_with_message('base_stack_1["base_ID1"]!=base_stack_2["base_ID1"]')

	if(base_stack_1["base_ID2"]!=base_stack_2["base_ID2"]): error_exit_with_message('base_stack_1["base_ID2"]!=base_stack_2["base_ID2"]')

	if(base_stack_1["face1"]!=base_stack_2["face1"]): return False

	if(base_stack_1["face2"]!=base_stack_2["face2"]): return False

	return True


########################################################################
def check_in_base_pair_list( template_base_pair, base_pair_list):

	for n in range(len(base_pair_list)):

		if(Is_equivalent_base_pair(template_base_pair, base_pair_list[n]) ): return True

	return False

########################################################################
def check_in_base_stack_list( template_base_stack, base_stack_list ):

	for n in range(len(base_stack_list)):

		if(Is_equivalent_base_stack(template_base_stack, base_stack_list[n]) ): return True

	return False

########################################################################

def Is_equivalent_base_pair(base_pair_1, base_pair_2, assert_valid=True):

	if(assert_valid):
		assert_valid_Base_pair(base_pair_1)
		assert_valid_Base_pair(base_pair_2)

	return (Is_equivalent_base_pair_strict(base_pair_1, base_pair_2) ) #Does not treat flip basepair as equivalent!

########################################################################

def Is_equivalent_base_pair_strict(base_pair_1, base_pair_2):

	if(base_pair_1["res1"]!=base_pair_2["res1"]): return False

	if(base_pair_1["res2"]!=base_pair_2["res2"]): return False

	if(base_pair_1["base_ID1"]!=base_pair_2["base_ID1"]): error_exit_with_message('base_pair_1["base_ID1"]!=base_pair_2["base_ID1"]')

	if(base_pair_1["base_ID2"]!=base_pair_2["base_ID2"]): error_exit_with_message('base_pair_1["base_ID2"]!=base_pair_2["base_ID2"]')

	if(base_pair_1["edge1"]!=base_pair_2["edge1"]): return False

	if(base_pair_1["edge2"]!=base_pair_2["edge2"]): return False

	if(base_pair_1["orientation"]!=base_pair_2["orientation"]): return False

	return True

########################################################################
def Is_canonical_BP(base_pair): #Watson-Crick or GU wobble!


	if(base_pair["edge1"] != "W"): return False
	if(base_pair["edge2"] != "W"): return False
	if(base_pair["orientation"] != "c"): return False

	base_ID1=base_pair["base_ID1"]
	base_ID2=base_pair["base_ID2"]

	assert_valid_base_ID(base_ID1)
	assert_valid_base_ID(base_ID2)

	if(base_ID1=="G" and base_ID2=="C"): return True
	if(base_ID1=="C" and base_ID2=="G"): return True
	if(base_ID1=="G" and base_ID2=="U"): return True
	if(base_ID1=="U" and base_ID2=="G"): return True
	if(base_ID1=="A" and base_ID2=="U"): return True
	if(base_ID1=="U" and base_ID2=="A"): return True

	return False

def Is_GU_wobble(base_pair):

	if(Is_canonical_BP(base_pair)==False): return False


	base_ID1=base_pair["base_ID1"]
	base_ID2=base_pair["base_ID2"]

	if(base_ID1=="G" and base_ID2=="U"): return True
	if(base_ID1=="U" and base_ID2=="G"): return True


	return False


########################################################################

def get_base_pair_recovery_stat( native_base_pair_list, model_base_pair_list , virtual_res_list, Verbose=False): #virtual_res are ignored!

	num_native_WC=0; 
	num_native_NWC=0;
	recovered_WC=0;
	recovered_NWC=0;

	for n in range(len(native_base_pair_list)):

		native_BP=native_base_pair_list[n]

		if(native_BP["res1"] in virtual_res_list): 
			print "Problem base_pair=", native_BP, " | virtual_res_list= ", virtual_res_list
			error_exit_with_message('native_BP["res1"] in virtual_res_list')

		if(native_BP["res2"] in virtual_res_list): 
			print "Problem base_pair=", native_BP, " | virtual_res_list= ", virtual_res_list
			error_exit_with_message('native_BP["res2"] in virtual_res_list')

		if( Is_canonical_BP(native_BP ) ):

			num_native_WC+=1						

			if ( check_in_base_pair_list( native_BP , model_base_pair_list ) ): recovered_WC+=1

		else:

			num_native_NWC+=1;


			if ( check_in_base_pair_list( native_BP, model_base_pair_list ) ):
				recovered_NWC+=1;
			else:
				if(Verbose): print "Missing native base_pair ", native_BP, " !"

	return (num_native_WC, num_native_NWC, recovered_WC, recovered_NWC)

########################################################################

def get_base_stack_recovery_stat(native_base_stack_list, model_base_stack_list , virtual_res_list, Verbose=False):

	num_native_BS=0; 
	recovered_BS=0;

	for n in range(len(native_base_stack_list)):

		native_BS=native_base_stack_list[n]

		if(native_BS["res1"] in virtual_res_list): 
			"Problem base_stack=", native_BS, " | virtual_res_list= ", virtual_res_list
			error_exit_with_message('native_BS["res1"] in virtual_res_list')

		if(native_BS["res2"] in virtual_res_list): 
			"Problem base_stack=", native_BS, " | virtual_res_list= ", virtual_res_list
			error_exit_with_message('native_BS["res2"] in virtual_res_list')

		num_native_BS+=1

		if( check_in_base_stack_list( native_BS, model_base_stack_list ) ):
			recovered_BS+=1;
		else:
			if(Verbose): print "Missing native base_stack ", native_BS, " !"

	return (num_native_BS, recovered_BS)


