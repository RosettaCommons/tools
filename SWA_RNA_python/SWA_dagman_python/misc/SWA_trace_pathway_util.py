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

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options



def get_parent_tag_and_file(desired_tag, infile, round_per_step):

	score_line_1=""
	score_line_2=""
	intermediate_silent_file=""
	score_lines=[]


	print "------------------------------------------------------------------------------------------------------"
	print  "desired_tag= %s // infile= %s " %(desired_tag, infile)

	for round_ID in range(round_per_step):

		num_tag_found=0
		found_parent_tag=False
		found_parent_file=False

		parent_tag=""
		parent_file=""

		if(exists(infile)==False): error_exit_with_message("infile(%s) doesn't exist" %(infile) )

		line_number=0

		for line in open(infile,'r'):

			line_number+=1

			if(round_ID==(round_per_step-1)): #last round
				if(line_number<=2): #the first two lines of the file
					score_lines.append(line)					

			if(num_tag_found>1): error_exit_with_message('num_tag_found>1') 

			if(num_tag_found==1 and found_parent_tag==False): #REMARK PARENT_TAG S_2_1
				line_list=line.split()
				if(line_list[0]!="REMARK"): error_exit_with_message('line_list[0]!="REMARK", line_list=%s' %(list_to_string(line_list) ) )
				if(line_list[1]!="PARENT_TAG"): error_exit_with_message('line_list[1]!="PARENT_TAG, line_list=%"' %(list_to_string(line_list) ) ) 
				found_parent_tag=True
				parent_tag=line_list[2]
				continue

			if(num_tag_found==1 and found_parent_file==False and found_parent_tag==True):  #REMARK SOURCE region_2_1_sample.cluster.out
				line_list=line.split()
				if(line_list[0]!="REMARK"): error_exit_with_message('line_list[0]!="REMARK", line_list=%s' %(list_to_string(line_list) ) )
				if(line_list[1]!="SOURCE"): error_exit_with_message('line_list[1]!="SOURCE, line_list=%"' %(list_to_string(line_list) ) ) 
				found_parent_file=True
				parent_file=line_list[2]
				continue


			if(line.find('SCORE:') == -1): continue  #Line doesn't contain the word SCORE
		
			if(line.find('desc') != -1): continue #Line is column_name description line
		
			tag=line.split()[-1] #Oct 19, 2010....assume that description is located at last column..change to this because some silent_file have structs with different score column length

			if(tag==desired_tag):
				if(round_ID==0):
					score_line_1=line
				else:
					score_line_2=line
			
				if(round_ID==(round_per_step-1)): score_lines.append(line)
			
				num_tag_found+=1

		if(num_tag_found==0):error_exit_with_message('could not find desired_tag (%s) in infile (%s)' %(desired_tag, infile) ) 

		#This only effect the final step:
		#the infile is the clustered silent_file. This converts S_2_1 to S_2 (the _1 is a modification done internally by the clusterer.
		if(round_per_step==2 and round_ID==0): 
			parent_tag_list=parent_tag.split("_")
			parent_tag=parent_tag_list[0] + "_" + parent_tag_list[1]

		print "round_ID= %s // parent_tag= %s // parent_file= %s" %(round_ID, parent_tag, parent_file)

		if(round_ID==1): intermediate_silent_file=infile

		desired_tag=parent_tag
		infile=parent_file

	if(round_per_step==2): #only for the FINAL step .. additional clustering step to create region_FINAL.out
		if(list_to_string(score_line_1.split()[:-1])!=list_to_string(score_line_2.split()[:-1])):
			print "list_to_string(score_line_1.split()[:-1])!=list_to_string(score_line_2.split()[:-1])"
			print "score_line_1: %s" %(score_line_1)
			print "score_line_2: %s" %(score_line_2)
			error_exit_with_message("list_to_string(score_line_1.split()[:-1])!=list_to_string(score_line_2.split()[:-1])")
	
	return (parent_tag, parent_file, score_lines, intermediate_silent_file)

