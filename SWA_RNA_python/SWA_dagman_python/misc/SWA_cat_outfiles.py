#!/usr/bin/env python

from sys import argv,exit,stdout,stderr
import sys
from os.path import exists
import string
import time
######################################################################

from SWA_dagman_python.utility.SWA_util import *
######################################################################

######################################################################
def concatenate_outfiles(infile_list, outfile, add_file_num_to_tag=True):

	print "outfile=%s | add_file_num_to_tag=%s " %(outfile, add_file_num_to_tag)
	fid = open( outfile, "w" )

	if(len(infile_list)==0):
		error_exit_with_message("len(infile)==0")

	start_concatenating_time=time.time()
##################(SPECIAL CASE: FIRST FILE UNIQUE OUTPUT)###########################################################

	if(exists(infile_list[ 0 ])==False): error_exit_with_message("infile_list[ 0 ] (%s) doesn't exist!")

	if(Is_valid_non_empty_silent_file(silent_file=infile_list[ 0 ])==False):  
		error_exit_with_message("infile_list[ 0 ] (%s) is not a valid non_empty_silent_file!" %(infile_list[ 0 ] ) )

	try:
		data = open( infile_list[ 0 ], 'r')
	except:
		error_exit_with_message("cannot open %s " %infile_list[0] )

	print 'writing file       #1: %s' %(infile_list[ 0 ])


	SEQUENCE_LINE = data.readline() # The SEQUENCE: gggcgcagccu line
	COLUMN_NAME_LINE   = data.readline() # The column name line
	COL_NAME_LIST=COLUMN_NAME_LINE.split()

	if(SEQUENCE_LINE[0:9]!='SEQUENCE:'):  		 error_exit_with_message("SEQUENCE_LINE[0:9]!='SEQUENCE:' for SEQUENCE_LINE (%s)" %(SEQUENCE_LINE) )
	if(COLUMN_NAME_LINE[0:6]!='SCORE:'):  		 error_exit_with_message("COLUMN_NAME_LINE[0:6]!='SCORE:'" )
	if(COLUMN_NAME_LINE.find('description') == -1 ):  error_exit_with_message("COLUMN_NAME_LINE.find('description') == -1" )
	if(COL_NAME_LIST[0]!='SCORE:'): error_exit_with_message("COL_NAME_LIST[0]!='SCORE:" )
	if(COL_NAME_LIST[-1]!='description'): error_exit_with_message("COL_NAME_LIST[-1]!='description'" )

	fid.write(SEQUENCE_LINE)
	fid.write(COLUMN_NAME_LINE)

	for line in data:
		fid.write( line )

	data.close()

####################ANOTHER FILE#####################################################################
	for i in range(1, len(infile_list)):

		if(i % 20 == 0):
			time_taken=time.time()-start_concatenating_time
			print "concatenated %d files in %f seconds" %(i, time_taken) 
			sys.stdout.flush()
			sys.stderr.flush()	

		if(exists( infile_list[i] )==False): error_exit_with_message("infile_list[%d] (%s) doesn't exist! " %(i, infile_list[i]) )
	
		if(Is_valid_non_empty_silent_file( silent_file=infile_list[ i ])==False):  
			error_exit_with_message("infile_list[ %d ] (%s) is not a valid non_empty_silent_file!" %(i, infile_list[ i ] ) )
		


		try:
			data = open(infile_list[i],'r')
		except:
			error_exit_with_message("cannot open %s " %(infile_list[i]) )
		
		print 'concatenating file #%d: %s' %(i+1, infile_list[i])

		first_line = data.readline() 
		second_line = data.readline()

		if(SEQUENCE_LINE!=first_line): 
			print "file_name=%s" %(infile_list[i] )
			print "SEQUENCE_LINE =", SEQUENCE_LINE
			print "first_line    =", first_line
			error_exit_with_message("SEQUENCE_LINE !=first_line")
 
		if( Is_equivalent_list(COL_NAME_LIST, second_line.split())==False ) : 
			print "file_name=%s" %(infile_list[i] )
			print "COL_NAME_LIST      =", COL_NAME_LIST
			print "second_line.split()=", second_line.split()
			error_exit_with_message("Is_equivalent_list(COL_NAME_LIST, second_line.split()==False)" )


		previous_line="INIT_PREVIOUS_LINE"


		while(True):
			line = data.readline()

			if(line==''): break #End of file!

			line = line[:-1]  ##remove the newline /n

			if( line.find( "REMARK BINARY_SILENTFILE" ) != -1 ):   
				if(line[0:24]!="REMARK BINARY_SILENTFILE"): error_exit_with_message("line[0:24]!=\"REMARK BINARY_SILENTFILE\"")
				continue	 #IGNORE REMARK_BINARY_SILENTFILE line.

			if(len(line) < 1): 
				print "filename=%s, previous_line= %s" %(infile_list[i], previous_line)
				error_exit_with_message("len(line) < 1!")
			previous_line=line

			#At this point, the line should be a data line or a REMARK line (REMARK PARENT_TAG, REMARK SOURCE)
			if(line.find( "REMARK ") != -1): #A REMARK LINE

				if(line[0:7]!="REMARK "): error_exit_with_message('line.find( "REMARK ") != -1 but line[0:7]!="REMARK ")| line=%s' %(line))

			else:

				#data line should contains a tag at the END of the line
				#tag should contains the "_" character. 
				cols = string.split( line )
				tag=cols[-1]
				if(tag.find("_") == -1): error_exit_with_message('tag.find("_") == -1 for tag (%s)' %(tag) )	

				###Move this inside the IF STATEMENT ON OCT 24, 2011.SHOULD NOT RENAME REMARK LINES!
				###assume that description (i.e. tag) is located at the last index
				if(add_file_num_to_tag):
					description_index = line.rfind(' ')+1 
					tag = line[description_index:]	
					newtag= '%d_%s' % (i, tag)		
					line = line[:description_index] + newtag
				################################################################################

			fid.write( line+'\n' )


		data.close()

	fid.close()

	if(add_file_num_to_tag==False): ensure_no_duplicate_tags(outfile)
