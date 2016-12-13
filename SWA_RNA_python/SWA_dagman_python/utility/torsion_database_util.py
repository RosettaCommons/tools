#!/usr/bin/env python

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
from os import system
from os.path import basename, dirname, exists, expanduser, abspath
from time import sleep
import popen2
import copy

######################################################################
from SWA_util import *

######################################################################


##########################################################################################################################################################

def extract_native_torsions(dope_in_native_torsions, native_torsions_only, START_torsion_database, native_pdb, All_sample_res_list):

	print "------------------ENTER extract_native_torsions-------------------------"
	print "dope_in_native_torsions=%s" %(dope_in_native_torsions)
	print "native_torsions_only=%s" %(native_torsions_only)
	print "START_torsion_database=%s" %(START_torsion_database)
	print "native_pdb=%s" %(native_pdb)
	print "All_sample_res_list=%s" %(All_sample_res_list)

	if(dope_in_native_torsions==True and native_torsions_only==True): error_exit_with_message("dope_in_native_torsions and native_torsions_only cannot both be True!")

	if(dope_in_native_torsions): 
		NEW_torsion_database="native_doped_" + START_torsion_database
	elif(native_torsions_only):	
		NEW_torsion_database="native_only.torsions" 	
	else:
		error_exit_with_message("Both dope_in_native_torsions and native_torsions_only equals False!")

	print "NEW_torsion_database=%s" %(NEW_torsion_database)

	if(exists(NEW_torsion_database)): error_exit_with_message("NEW_torsion_database (%s) already exist!" %(NEW_torsion_database) )

	TORSION_DB=open( NEW_torsion_database, 'w')

	#create the native_torsion file.
	submit_subprocess("extract_torsions_from_pdb.py -pdb %s -sample_res %s " %(native_pdb, list_to_string(All_sample_res_list) ) )
	
	#check that the output torsion files exist.
	native_torsion_file="sample_res_" + native_pdb.replace(".pdb", "") + ".torsions"
	idealized_native_torsion_file="sample_res_" + "idealize_" + native_pdb.replace(".pdb", "") + ".torsions"

	if(exists(native_torsion_file)==False): error_exit_with_message("native_torsion_file (%s) doesn't exist!" %(native_torsion_ofile) )
	if(exists(idealized_native_torsion_file)==False): error_exit_with_message("idealized_native_torsion_file (%s) doesn't exist!" %(idealized_native_torsion_file) )
		
	
	#could get use cat command...but this is probably safer.
	infile_torsion=open( native_torsion_file ).readlines()
	for line in infile_torsion: TORSION_DB.write(line)

	infile_torsion=open( idealized_native_torsion_file ).readlines()
	for line in infile_torsion: TORSION_DB.write(line)

	if(dope_in_native_torsions):
		infile_torsion=open( START_torsion_database ).readlines()
		for line in infile_torsion: TORSION_DB.write(line)

	TORSION_DB.close()

	print "------------------EXIT extract_native_torsions-------------------------" 

	return NEW_torsion_database

#########################################################################################################################################

def create_excised_torsion_database(torsion_database, excise_segment_torsion_DB, ensure_excise_segment_seq_match, sequence, sample_res_list):

	print "----------Enter create_excised_torsion_database-----------"	

	print "excise_segment_torsion_DB=%s" %(excise_segment_torsion_DB) 
	print "ensure_excise_segment_seq_match=%s" %(ensure_excise_segment_seq_match)

	TORSION_DB=open( torsion_database).readlines()

	excised_torsion_database="excised_segment" 
	
	excised_line_pair_list=[]

	for excise_segment_string in excise_segment_torsion_DB:

		print "excise_segment_string=%s" %(excise_segment_string) 

		excise_segment=excise_segment_string.split("-")

		for n in range(len(excise_segment)):
			excise_segment[n]=int(excise_segment[n])
		
		if(len(excise_segment)!=2):
			error_exit_with_message("len(excise_segment)!=2, excise_segment= %s " %( list_to_string(excise_segment) ) )

		if(excise_segment[0]>=excise_segment[1]):
			error_exit_with_message("excise_segment_torsion_DB[0]=%d=>%d=excise_segment[1] " %(excise_segment[0], excise_segment[1]) )

		excised_torsion_database+="_%d_%d" %( excise_segment[0], excise_segment[1] )

		excise_segment_length=(excise_segment[1]-excise_segment[0])+1

		#############################consistency check##########################################################
		sampling_base_ID=[]

		if(ensure_excise_segment_seq_match): #This assumes that sample_res_list is a single continuous region and that there is only one excise segment.

			sample_res_list.sort()

			for n in range(len(sample_res_list)):
				if(n==0):continue
				if( (sample_res_list[n])!=(sample_res_list[n-1]+1) ):
					error_exit_with_message("(sample_res_list[n])=(%d)!=(%d)=(sample_res_list[n-1]+1)" %(sample_res_list[n], sample_res_list[n-1]+1) )

			for n in range(len(sample_res_list)+2):

				seq_num= n + (sample_res_list[0]-1)

				#print "native_pdb_seq_num= %d " %(seq_num)

				sampling_base_ID.append(sequence[seq_num-1])

			print "sampling_base_ID= ", sampling_base_ID
	
			if(len(sampling_base_ID)!=excise_segment_length):
				error_exit_with_message("len(sampling_base_ID)=(%d)!=(%d)=excise_segment_length" %( sampling_base_ID, excise_segment_length ) )
		#######################################################################################
			
		start_excised_line=0
		end_excised_line=0

		num_excised_line_region=0

		for line_num in range(len(TORSION_DB)):

			start_seq_num=int(TORSION_DB[line_num].split()[-1])

	#		print "TEST seq_num= %d " %( start_seq_num ) 
		
			if(start_seq_num==excise_segment[0]):

				print "potential_excise_segment match starting at line_num=%d" %(line_num)

				increment=0
				while(True):
			
					seq_num=int(TORSION_DB[line_num+increment].split()[-1])										
					base_type=TORSION_DB[line_num+increment].split()[0]	
		
					print "increment=%d , line_num=%d, seq_num=%d , base_type=%s " %( increment, line_num+increment, seq_num, base_type  ) 

					if(seq_num!=start_seq_num+increment):
						print "seq_num!=start_seq_num+increment, seq_num=%d, start_seq_num=%d, increment=%d " %( seq_num, start_seq_num, increment )
						break

					if(ensure_excise_segment_seq_match):
						if(base_type!=sampling_base_ID[increment]):
							print "base_type!=sampling_base_ID[increment], base_type=%s, sampling_base_ID[increment]=%s, increment=%d " %(  base_type, sampling_base_ID[increment], increment)
							break
	

					if(increment>=excise_segment_length):
						error_exit_with_message("increment=(%d)>(%d)=excise_segment_length" %(increment, excise_segment_length) )

					if( increment==(excise_segment_length-1) ):	
						num_excised_line_region+=1
						start_excised_line=line_num
						end_excised_line=line_num+increment
						print "FOUND excise_segment match, num_excised_line_region_so_far=%d, start_excised_line=%d, end_excised_line=%d" %(num_excised_line_region, start_excised_line, end_excised_line)
						break

					increment+=1

		if(num_excised_line_region!=1): error_exit_with_message("num_excised_line_region (%d) != 1" %(num_excised_line_region) )

		excised_line_pair_list.append([start_excised_line, end_excised_line])

	#########################################################

	excised_torsion_database+="_%s" %(basename(torsion_database))

	print "excised torsion_database= %s " %(excised_torsion_database)			

	EXCISED_TORSION_DB=open( excised_torsion_database,  'w')

	for line_num in range( len(TORSION_DB) ):

		seq_num=int(TORSION_DB[line_num].split()[-1])										
		base_type=TORSION_DB[line_num].split()[0]	
		
		Is_excise_line=False	

		for excised_line_pair in excised_line_pair_list:
			if( line_num>=excised_line_pair[0] and line_num<=excised_line_pair[1] ):
				print "excising_line_num=%d, seq_num=%d, base_type=%s " %(line_num,seq_num, base_type)
				Is_excise_line=True

		if(Is_excise_line): continue

		EXCISED_TORSION_DB.write( TORSION_DB[line_num] )
		

	EXCISED_TORSION_DB.close()




	print "----------Exit create_excised_torsion_database-----------"	

	return excised_torsion_database

