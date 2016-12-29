#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
import os
import copy
from error_util import *






####################################################################

def list_to_string(my_list, separator=" ", first_separator=True):

	#my_list=my_list #necessary? what does did do? Mod out on July 31, 2011.

	if(isinstance( my_list, list )==False):
		print "problem my_list:", my_list
		error_exit_with_message("my_list object is not a list!")

	my_string=""

	for n in range(len(my_list)):

		if(n==0 and (first_separator==False) ):
			my_string+=str(my_list[n])
		else:
			my_string+=separator + str(my_list[n])

	return my_string

####################################################################
def assert_no_duplicate_in_string_list(string_list): #To should work for both python and Rosetta (C++) commands

	for str_element in string_list:

		if(isinstance( str_element , str )==False):
			print "\n\nERROR string_list=", string_list
			error_exit_with_message("string_list is not a str!")

	str_elements_so_far=[]

	for n in range(len(string_list)):

		str_element=string_list[n]

		if(str_elements_so_far.count(str_element)>0):
			exit_message ="\n\nERROR string_list=%s \n" %(list_to_string(string_list))
			exit_message+="str_element (%s) already exist in the str_elements_so_far ( %s )\n" %(str_element, list_to_string(str_elements_so_far, " | ") )
			error_exit_with_message(exit_message)

		str_elements_so_far.append(str_element)

################################################################################################################
def add_boundary_seq_num_to_list(start_list, max_seq_num, boundary_size):

	list_with_boundary=[]

	if(isinstance( start_list, list )==False): error_exit_with_message("start_list object is not a list!")
	if(isinstance( max_seq_num, int )==False): error_exit_with_message("max_seq_num object is not an int!")
	if(isinstance( boundary_size, int )==False): error_exit_with_message("boundary_size object is not an int!")

	if(len(start_list)==0): error_exit_with_message("len(start_list)==0!")

	for start_seq_num in start_list:
		if(isinstance( start_seq_num, int )==False): error_exit_with_message("start_seq_num object is not an int!")

		if(start_seq_num<1): error_exit_with_message("start_seq_num (%s) is less than zero!" %(start_seq_num))

		if(start_seq_num>max_seq_num): error_exit_with_message("start_seq_num (%s) greater than max_seq_num (%s)" %(start_seq_num, max_seq_num))

	for seq_num in range(1, max_seq_num+1):

		if(seq_num in list_with_boundary): continue

		include_seq_num=False

		#FIX ON March 12, 2012: USED TO BE if((seq_num) in list_with_boundary): include_seq_num=True
		#THE PREVIOUS CASE LEAD TO A BUG IF THERE IS A NON-CONSECUTIVE SEQ_NUM!
		#HENCE ONLY synG_antiG_bp and short_Human_HAR1, short_Chimp_HAR1 ARE EFFECT [mainly in DAG_rebuild_bulge.py + calc_chem_shift]

		if(seq_num in start_list): include_seq_num=True

		for n in range(1, boundary_size+1):

			if( (seq_num-n) in start_list): include_seq_num=True

			if( (seq_num+n) in start_list): include_seq_num=True

		if(include_seq_num): list_with_boundary.append(seq_num)

	list_with_boundary.sort()

	return list_with_boundary

################################################################################################################


def list_intersection(list_one, list_two):

	if(isinstance( list_one, list )==False): error_exit_with_message("list_one object is not a list!")
	if(isinstance( list_two, list )==False): error_exit_with_message("list_two object is not a list!")

	intersection=[]

	for seq_num_one in list_one:
		if(seq_num_one in list_two): intersection.append(seq_num_one)

	intersection=list(Set(intersection))

	intersection.sort()

	return intersection

################################################################################################################

def Is_disjoint_list(list_one, list_two):

	if(isinstance( list_one, list )==False): error_exit_with_message("list_one object is not a list!")
	if(isinstance( list_two, list )==False): error_exit_with_message("list_two object is not a list!")

	return (len(list_intersection(list_one, list_two) )==0 )

################################################################################################################


def Is_subset_list(smaller_list, bigger_list):

	if(isinstance( smaller_list, list )==False): error_exit_with_message("smaller_list object is not a list!")
	if(isinstance( bigger_list, list )==False): error_exit_with_message("bigger_list object is not a list!")

	if( len(list( Set(smaller_list)-Set(bigger_list) ) )==0 ):
		return True
	else:
		return False

############################################################################################################

def get_segment_string_list(res_list, cutpoint_open_list=[]):

	if(isinstance( res_list, list )==False):
		print "PROBLEM res_list:", res_list
		error_exit_with_message("res_list object is not a list!")

	for seq_num in res_list:
		if(isinstance( seq_num, int )==False):
			print "PROBLEM res_list:", res_list
			error_exit_with_message("seq_num object is not an int!")

	segment_list=[]
	found_start=False
	found_end=False

	for n in range(len(res_list)):

		curr_seq_num=res_list[n]

		if(found_start==False):
			start_res=curr_seq_num
			found_start=True

		if(found_start==True):
			if( (n== (len(res_list)-1)) or ((curr_seq_num+1)!=res_list[n+1]) or (curr_seq_num in cutpoint_open_list) ):
				found_end=True



		if(found_end==True):
			if (start_res == curr_seq_num):
				segment_list.append("%d" %(curr_seq_num) )
			else:
				segment_list.append("%d-%d" %(start_res, curr_seq_num) )
			found_start=False
			found_end=False

	return segment_list

############################################################################################################
def seq_num_list_to_string(res_list):

	if(isinstance( res_list, list )==False):
		print "PROBLEM res_list:", res_list
		error_exit_with_message("res_list object is not a list!")

	for seq_num in res_list:
		if(isinstance( seq_num, int )==False):
			print "PROBLEM res_list:", res_list
			error_exit_with_message("seq_num object is not an int!")

	segment_list=get_segment_string_list(res_list)

	seq_num_string=""

	for segment in segment_list:

		if(len(segment.split("-"))!=2): error_exit_with_message("len(segment.split(\"-\"))!=2")

		start_res=int(segment.split("-")[0])
		end_res=  int(segment.split("-")[1])

		if(start_res==end_res):
			seq_num_string+="_%s" %(start_res)
		else:
			seq_num_string+="_%s_to_%s" %(start_res, end_res)

	return seq_num_string

############################################################################################################
def parse_segment_string_list(segment_string_list, verbose=False, do_sort=True):

	if(verbose): print "parsing seg_num_list from segment_string_list=" , segment_string_list
	if(segment_string_list==['']): return []

	seq_num_list=[]

	for segment_string in segment_string_list:

		if(len(segment_string.split("-"))==1): #single seq_num
			seq_num_list.append(int(segment_string))
		elif(len(segment_string.split("-"))==2):
			start_segment=int(segment_string.split("-")[0])
			end_segment=int(segment_string.split("-")[1])
			if(start_segment>end_segment): error_exit_with_message("start_segment(%d)>end_segment(%d)" %(start_segment, end_segment))
			seq_num_list.extend(range(start_segment, end_segment+1))
		else:
			print "segment_string_list=", segment_string_list
			error_exit_with_message("ERROR segment_string=%s, len(segment_string.split(\"-\"))=%d" %(segment_string, len( segment_string.split("-") ) ) )
		if(do_sort): seq_num_list.sort()

	if( len(seq_num_list)!=len(list(Set(seq_num_list) ) ) ):
		print "ERROR:      segment_string_list=", segment_string_list
		print "ERROR:             seq_num_list=", seq_num_list
		print "ERROR: list(Set(seq_num_list) )=",  list(Set(seq_num_list) )
		error_exit_with_message("len(seq_num_list)(%d)!=(%d)len(list(Set(seq_num_list) ) )" %( len(seq_num_list), len(list(Set(seq_num_list) ) ) ) )

	return seq_num_list

############################################################################################################
def Is_equivalent_list(list_one, list_two):

	if( len(list_one)!=len(list_two) ): return False

	for n in range(len(list_one)):
		if(list_one[n]!=list_two[n]): return False


	return True


####################################################################

def print_seq_num_list( INPUT_list_text_name, INPUT_seq_num_list, INPUT_text_spacing=30):

	if(isinstance( INPUT_list_text_name, str )==False): error_exit_with_message("INPUT_list_text_name object is not an str!")

	if(isinstance( INPUT_seq_num_list, list )==False): error_exit_with_message("INPUT_list_text_name object is not an int!")

	for seq_num_obj in INPUT_seq_num_list:
		if( isinstance( seq_num_obj, int)==False):  error_exit_with_message("seq_num_obj object is not an int!")

	if(isinstance( INPUT_text_spacing, list )==False): error_exit_with_message("INPUT_text_spacing object is not an list!")


	text_spacing=copy.deepcopy(INPUT_text_spacing)

	list_text_name=copy.deepcopy(INPUT_list_text_name)

	seq_num_list=copy.deepcopy(INPUT_seq_num_list)

	seq_num_list.sort()

	#Variable string length in PYTHON:
	#("%i" % 12).rjust(a)
	#Or, more ugly:
	#"%%%di" % a % 12

	screen_width=154
	seq_num_width=3

	text_width=max(INPUT_text_spacing, len(INPUT_list_text_name) )

	curr_string_length=text_width

	print ("%s" % INPUT_list_text_name).rjust(INPUT_text_spacing),

	seq_num=1;
	for n in range(len(seq_num_list)):

		while(seq_num<seq_num_list[n]):
			print ("%s" %("")).rjust(seq_num_width),
			curr_string_length+=seq_num_width
			if(curr_string_length>(screen_width-seq_num_width-4)):
				curr_string_length=0
				print "cont."
				print ("").rjust(text_width),

			seq_num+=1

		print ("%s" %(seq_num_list[n])).rjust(seq_num_width),
		seq_num+=1;
		curr_string_length+=seq_num_width
		if(curr_string_length>(screen_width-seq_num_width-4)):
			curr_string_length=0
			print "cont."
			print ("").rjust(text_width),


	print ""


####################################################################
def print_condense_seq_num_list(prestring, seq_num_list):

	if(isinstance( prestring, str )==False): error_exit_with_message("prestring object is not a str!")

	if(isinstance( seq_num_list, list )==False): error_exit_with_message("seq_num_list object is not an list!")

	for seq_num in seq_num_list:
		if( isinstance( seq_num, int)==False):  error_exit_with_message("seq_num object is not an int!")

	print "%s %s" %(prestring, list_to_string( get_segment_string_list(seq_num_list) ) )

####################################################################
def list_to_string_with_dashes(seq_num_list):

	return string.join( get_segment_string_list( seq_num_list ) )


####################################################################

