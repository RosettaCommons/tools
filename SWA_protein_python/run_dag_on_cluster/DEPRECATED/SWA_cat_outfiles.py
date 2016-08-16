#!/usr/bin/env python

from sys import argv,exit,stdout,stderr
import sys
from os.path import exists
import string
from SWA_util import *


def concatenate_outfiles(infile_list, outfile):

	print "outfile= ", outfile
	fid = open( outfile, "w" )

	if(len(infile_list)==0):
		error_exit_with_message("len(infile)==0")

##################(SPECIAL CASE: FIRST FILE UNIQUE OUTPUT)###########################################################

	try:
		data = open( infile_list[ 0 ], 'r')
	except:
		error_exit_with_message("cannot open %s " %infile_list[0] )

	for line in data:
		fid.write( line )

	data.close()

####################ANOTHER FILE#####################################################################
	for i in range(1, len(infile_list)):

		if not exists( infile_list[i] ):
			print('Does not exist! '+infile_list[i] )
			sys.exit(1)
	
		try:
			data = open(infile_list[i],'r')
		except:
			error_exit_with_message("cannot open %s " %(infile_list[i]) )
		
		print 'concatenating %s' %(infile_list[i])
		line = data.readline() # The SEQUENCE: gggcgcagccu line
		line = data.readline() # The column name line

		while line:
			line = data.readline()[:-1]  ##remove the newline /n

			if (line.find( "REMARK BINARY_SILENTFILE" ) != -1): continue	 #This include  REMARK_BINARY_SILENTFILE
			if len(line) < 1: continue #Empty line

			###Consistency check...at this point, the line should be a data line or a REMARK line (REMARK PARENT_TAG, REMARK SOURCE)
			###A characteristic of the data line is that it contains a tag at the line of the line
			###A characteristic of the tag is that it contains the "_" character. 
			consistency_check=False

			if(line.find( "REMARK ") != -1): #REMARK LINE
#				print "REMARK LINE: ", line
				consistency_check=True
			else: #find the tag
				cols = string.split( line )
				tag=cols[-1]
#				print "tag: ", tag
				if(tag.find("_") != -1): consistency_check=True

			if(consistency_check==False):
				error_exit_with_message("Fail consistency check" )				

			#assume that description is located at the last index
			description_index = line.rfind(' ')+1 

			tag = line[description_index:]	
			newtag= '%d_%s' % (i, tag)		
			line = line[:description_index] + newtag
	#		print "newtag=" + newtag
			fid.write( line+'\n' )
		
		data.close()

	fid.close()
