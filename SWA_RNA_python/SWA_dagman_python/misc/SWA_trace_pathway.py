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

from SWA_trace_pathway_util import *
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################


#Apr_2_TRAIL_2_MANUAL_KEEP_LOG_FILES

#SWA_trace_pathway.py -tag S_2 -silent_file region_FINAL.out -num_elements 10 -rmsd_column 27 >  SWA_trace_pathway_S_2_region_FINAL.txt

#SWA_trace_pathway.py -tag S_5 -silent_file  region_0_9_sample.cluster.out -num_elements 10 -rmsd_column 27 >  SWA_trace_pathway_S_5_region_0_9_sample.cluster.txt


#SWA_trace_pathway.py -tag S_1 -silent_file region_1_0_sample.cluster.out -num_elements 6 -rmsd_column 27 >  SWA_trace_pathway_S_1_region_1_0_sample.cluster.txt

#SWA_trace_pathway.py -tag S_3 -silent_file region_1_0_sample.cluster.out -num_elements 6 -rmsd_column 27 -KEEP_LOG_FILE_FOLDER Apr_2_TRAIL_2_MANUAL_KEEP_LOG_FILES/ >  SWA_trace_pathway_S_3_region_1_0_sample.cluster.out.txt



#SWA_trace_pathway.py -tag  M_000017395757_5 -silent_file REGION_1_0/start_from_region_3_0_sample_filtered.out  -final_silent_is_not_clustered True -num_elements 6 -rmsd_column 27 -KEEP_LOG_FILE_FOLDER Apr_2_TRAIL_2_MANUAL_KEEP_LOG_FILES/ >  SWA_trace_pathway_S_0_REGION_1_0_SCORE_start_from_region_3_0_sample.out.txt


#SWA_trace_pathway.py -tag  S_0 -silent_file region_6_5_sample.cluster.out  -num_elements 10 -rmsd_column 27 >  SWA_trace_pathway_S_0_region_6_5_sample.cluster.txt


#SWA_trace_pathway.py -tag  S_77 -silent_file region_8_4_sample.cluster.out -num_elements 12 -stop_region_list 11-0 0-1 -rmsd_column 27 > S_77_region_8_4_path_log.txt

#SWA_trace_pathway.py -tag  S_1178 -silent_file region_FINAL.out -num_elements 12 -stop_region_list 11-0 0-1 -rmsd_column 26 > S_1178_path.txt

#SWA_trace_pathway.py -tag  S_0 -silent_file region_FINAL.out -num_elements 5 -rmsd_column 27 > S_0_path.txt

#SWA_trace_pathway.py -tag  S_6 -silent_file region_6_5_sample.cluster.out -num_elements 10  -rmsd_column 29 > S_6_region_6_5_path_log.txt


stop_region_list = parse_options( argv, "stop_region_list", [""] )
final_tag= parse_options( argv, "tag", "" )
final_silent_file= parse_options( argv, "silent_file", "" )
output_score_files= parse_options( argv, "output_score_files", "True" )
extract_pdb= parse_options( argv, "extract_pdb", "True" )
num_elements= parse_options( argv, "num_elements", 0)
rmsd_column=  parse_options( argv, "rmsd_column", 27 )

KEEP_LOG_FILE_FOLDER= parse_options( argv, "KEEP_LOG_FILE_FOLDER", "KEEP_LOG_FILE/")
final_silent_is_not_clustered= parse_options( argv, "final_silent_is_not_clustered", "False")

if(num_elements==0): error_exit_with_message('num_elements==0')  

if(final_tag==""): error_exit_with_message('final_tag==""')  
if(final_silent_file==""): error_exit_with_message('final_silent_file==""')  

print "------------------------------------------------------------------------------------------------------"
print "------------------------------------------------------------------------------------------------------"


desired_tag_list=[final_tag]
infile_list=[final_silent_file]
pdb_info_list=[]

step_ID=0

output_folder="TRACE_PATHWAY_final_%s_%s/" %(final_silent_file.replace("/", "__"), final_tag)

if(output_score_files):
	if(exists(output_folder)): submit_subprocess( "rm -r %s" %(output_folder) )
	submit_subprocess( "mkdir %s" %(output_folder) )

while(len(infile_list)>0):
	
	next_infile_list=[]
	next_desired_tag_list=[]

	if(len(infile_list)!=len(desired_tag_list)):  error_exit_with_message('len(infile_list)!=len(desired_tag_list)') 

	for file_num in range(len(infile_list)):

		infile=infile_list[file_num]
		desired_tag=desired_tag_list[file_num]

		print "------------------------------------------------------------------------------------------------------"
		print "------------------------------------------------------------------------------------------------------"
		print "STEP_ID= %d" %(step_ID)
		if(infile=="region_FINAL.out"):
			round_per_step=2
		else:
			round_per_step=1

		if(step_ID==0 and final_silent_is_not_clustered==True):
			parent_tag=desired_tag
			parent_file=infile
			score_lines=[]
			intermediate_silent_file=""
		else:
			(parent_tag, parent_file, score_lines, intermediate_silent_file)=get_parent_tag_and_file(desired_tag, infile, round_per_step)

			if(parent_file[0:3]=="../"): parent_file=parent_file[3:]

		if(output_score_files):
			output_filename=""
			if(step_ID==0 and final_silent_is_not_clustered==True):
				output_filename="%s/SCORE_%s" %( output_folder, infile.replace("/", "__") )
			elif(infile=="region_FINAL.out"):
				if(intermediate_silent_file==""): error_exit_with_message('infile=="region_FINAL.out" but intermediate_silent_file==""')
				output_filename="%s/SCORE_%s" %(output_folder,intermediate_silent_file)
			else:
				output_filename="%s/SCORE_%s" %(output_folder,infile)
	
			if(exists(output_filename)): error_exit_with_message("output_filename (%s) already exist!" %(output_filename) )
	
			print "output_filename= %s " %(output_filename)	

			OUTFILE=open( output_filename,'w')
			for score_line in score_lines:
				OUTFILE.write( "%s" %(score_line) )
			OUTFILE.close()

		pdb_info={}
		pdb_info["silent_file"]=infile
		pdb_info["tag"]=desired_tag
		pdb_info_list.append(pdb_info)

		if(infile=="region_FINAL.out"):
			inter_pdb_info={}
			inter_pdb_info["silent_file"]=intermediate_silent_file
			inter_pdb_info["tag"]=score_lines[-1].split()[-1]
			pdb_info_list.append(inter_pdb_info)
		##############################################################################

		queue_ID=0

		if(parent_tag[0]=='M'):
			queue_ID=0
		else:
			queue_ID=int(parent_tag.split('_')[0])	 #WARNING THIS is not actually the queue_ID and instead is just the access order of the file from the glob_string. 

		#print "queue_ID= %d" %(queue_ID)

		##############################################################################
		end_region=(parent_file.split('/')[0])
		start_region=( ( parent_file.split('/')[1]).replace('_sample_filtered.out','') ).upper() 

		# parent_file= REGION_3_2/start_from_region_3_1_sample_filtered.out
		# -->REGION_3_2_START_FROM_REGION_3_1

		reducer_log_file=KEEP_LOG_FILE_FOLDER + "/CONDOR/" +end_region + '_' + start_region + '/outfile/reducer.out'

		#print "reducer_log_file= %s " %(reducer_log_file)

		if(exists(reducer_log_file)==False): error_exit_with_message("reducer_log_file (%s) doesn't exist!" %(reducer_log_file))
		##############################################################################
		#i=27, globfiles[i]=REGION_3_2/START_FROM_REGION_3_1_S_82/region_3_2_sample.out

		num_queue_ID_found=0
		sampler_silent_file=""

		for line in open(reducer_log_file,'r'):
	
			if(line[0:2]!='i='): continue

			if(line.split()[1].split('=')[0] != 'globfiles[i]' ): continue
		
			if(queue_ID!=int((line.split()[0]).split('=')[1].replace(',','') ) ): continue
	
			num_queue_ID_found+=1
			sampler_silent_file=line.split()[1].split('=')[1]

		if(num_queue_ID_found!=1): error_exit_with_message("num_queue_ID_found (%s) !=1" %(num_queue_ID_found) ) 

		#print "sampler_silent_file= %s " %(sampler_silent_file)  #REGION_3_2/START_FROM_REGION_3_1_S_45/region_3_2_sample.out

		start_string_list=sampler_silent_file.split('/')[1].split('_')

		#print "start_string_list=" , start_string_list 


		combine_long_loop=False
		num_starting_infile=1

		if(len(start_string_list)>=6 and start_string_list[5]=="AND"):
			combine_long_loop=True
			num_starting_infile=2

		curr_silent_file_desired_tag_list=[]

		ACTUAL_queue_ID=int(start_string_list[-1])

		print "ACTUAL_queue_ID= %d" %(ACTUAL_queue_ID)

		if(combine_long_loop):			
			filter_out_file="FILTER_OUTFILE/%s_%s_combine_long_loop_filterer_output.txt" %(end_region.lower(), start_region.lower() )
			print "filter_out_file= %s" %(filter_out_file) 
			if(exists(filter_out_file)==False): error_exit_with_message("filter_out_file (%s) doesn't exist!" %(filter_out_file) )

			num_filterer_queue_ID_found=0
			filterer_info_line=""
			for line in open(filter_out_file,'r'):
				pose_num=int(line.split()[3])
				corr_queue_ID=pose_num-1
				if(corr_queue_ID==(ACTUAL_queue_ID)):
					filterer_info_line=line.split()
					num_filterer_queue_ID_found+=1
		
			if(num_filterer_queue_ID_found!=1): error_exit_with_message("num_filterer_queue_ID_found (%s) !=1" %(num_filterer_queue_ID_found) ) 
	
			print "filterer_info_line= %s " %(filterer_info_line)
			curr_silent_file_desired_tag_list.append(	filterer_info_line[0])		
			curr_silent_file_desired_tag_list.append(	filterer_info_line[1])	

			
			#/FILTER_OUTFILE/region_7_6_start_from_region_0_6_and_8_0_combine_long_loop_filterer_output.txt

		else:
			curr_silent_file_desired_tag_list.append("%s_%d" %(start_string_list[-2],ACTUAL_queue_ID) )

		for jj in range(num_starting_infile):

			if(jj==0): #append pose if combine_long_loop
				element_i=int(start_string_list[3])
				element_j=int(start_string_list[4])
			else: #prepend pose if combine_long_loop
				element_i=int(start_string_list[6])
				element_j=int(start_string_list[7])

			next_infile='%s_%s_%s_sample.cluster.out'  %(start_string_list[2].lower(), element_i, element_j)     # region_3_1_sample.cluster.out

			if(element_i==element_j):
				print "------------------------------------------------------------------------------------------------------"
				print "REACH STARTING POINT! element_i(%s)==element_j(%s)" %(element_i, element_j)
				print "next_infile %s" %(next_infile)
				continue

			###########################################################
			reach_stop_region=False
			if(stop_region_list!=[""]):
				for nn in range(len(stop_region_list)):
					stop_region=stop_region_list[nn].split('-')
					if(len(stop_region)!=2): error_exit_with_message("len(stop_region_list[%d])!=2, stop_region_list[%d]=" %(nn,nn,stop_region_list[nn]) )
					if( element_i==int(stop_region[0]) and element_j==int(stop_region[1]) ): 
						reach_stop_region=True
						break

			if(reach_stop_region==True): 
				print "------------------------------------------------------------------------------------------------------"
				print "REACH STOP REGION! element_i(%s)==element_j(%s)" %(element_i, element_j)
				print "next_infile %s" %(next_infile)
				continue
			###########################################################

			step_ID+=1
			print "next_desired_tag= %s // next_infile= %s " %(curr_silent_file_desired_tag_list[jj], next_infile)
			next_infile_list.append(next_infile)
			next_desired_tag_list.append(curr_silent_file_desired_tag_list[jj])
			
	infile_list=next_infile_list
	desired_tag_list=next_desired_tag_list
print "------------------------------------------------------------------------------------------------------"

if(output_score_files):
	gnuplot_command="create_gnuplot_script.py -trace_pathway_folder %s %s quick -rmsd_column %s " %(output_folder,num_elements, rmsd_column)
	print "gnuplot_command= %s " %(gnuplot_command) 
	submit_subprocess(gnuplot_command)

print "------------------------------------------------------------------------------------------------------"

base_folder=os.path.abspath(".")

os.chdir( output_folder )

for pdb_info in pdb_info_list:
	pdb_info["silent_file"]=os.path.abspath("../" + pdb_info["silent_file"])
	print pdb_info
	if(extract_pdb): 
		
		if(exists("%s.pdb" %(pdb_info["tag"]) ) ): error_exit_with_message('%s.pdb already exist!' %(pdb_info["tag"]) )
		submit_subprocess("SWA_extract_pdb.py -silent_file %s -tag %s -move_silent_file_into_folder False" %(pdb_info["silent_file"], pdb_info["tag"]) )
		submit_subprocess("mv %s.pdb %s_%s.pdb" %(pdb_info["tag"], pdb_info["tag"], basename(pdb_info["silent_file"]).replace(".out","") ) ) 


print "base_folder %s"  %(base_folder)
os.chdir( base_folder )

print "------------------------------------------------------------------------------------------------------"
print "------------------------------------------------------------------------------------------------------"
print "SWA_trace_pathway.py SUCCESSFULLY RAN!"



