#!/usr/bin/env python

######################################################################

from SWA_dagman_python.utility.SWA_util import *
######################################################################

from sys import argv,exit,stderr
from os import system
from os.path import basename
import string
import time
#import sets
from sets import Set

######################################################################



#print "max_score_gap=", max_score_gap
def filter_silent_file(silent_file, outfile_name, scorecol_name="score", REVERSE=True, max_n_struct=0, max_score_gap=0.0, abs_score_cut="FALSE", remove_SCORE_file=False):

	verbose=True

	title='filter_silent_file()'

	print_title_text("Enter " + title)

	print "silent_file=%s " %silent_file 
	print "outfile_name=%s " %outfile_name
	print "REVERSE= ", REVERSE 

	if(silent_file == ""): error_exit_with_message("need to input silent_file!")

	if(outfile_name == ""): error_exit_with_message("need to outfile_name")

	if(exists(outfile_name)): #Jan 10, 2012
		print "Warning outfile (%s) already exist.....removing!" %(outfile_name)
		submit_subprocess("rm %s" %(outfile_name))

	mode_max_struct=False
	mode_score_gap=False
	mode_abs_score_cut=False

	print "scorecol_name=%s" %(scorecol_name)

	num_mode=0

	if(max_n_struct !=0): 
		num_mode+=1
		mode_max_struct=True
		print "mode_max_struct=True! | max_n_struct=%d, max_score_gap=%8.3f" %(max_n_struct, max_score_gap)

	if(max_score_gap > 0.01): 
		num_mode+=1
		mode_score_gap=True
		max_n_struct=99999999999999999999
		if(max_score_gap<0.0): error_exit_with_message("max_score_gap= %8.3f need to be a positive real number" %(max_score_gap))
		print "mode_score_gap=True! | max_n_struct=%d, max_score_gap=%8.3f" %(max_n_struct, max_score_gap)

	if(abs_score_cut!="FALSE"): 
		num_mode+=1
		mode_abs_score_cut=True
		max_n_struct=99999999999999999999

		try:
			abs_score_cut=float(abs_score_cut)
		except:
			error_exit_with_message("Cannot convert abs_score_cut to float! abs_score_cut=%s" %(abs_score_cut))

		print "mode_abs_score_cut! | max_n_struct=%d, abs_score_cut=%8.3f" %(max_n_struct, abs_score_cut)

	if(num_mode!=1): error_exit_with_message("num_mode!=1")

	#######################################################################################################################################
	##Update on Jan 10, 2012
	try:
		silent_data = open( silent_file , 'r')
	except:
		error_exit_with_message("cannot open (%s) " %(silent_data))


	SEQUENCE_LINE 		= silent_data.readline() # The SEQUENCE: gggcgcagccu line
	COLUMN_NAME_LINE = silent_data.readline() # The column name line
	THIRD_LINE				=	silent_data.readline() #Possible REMARK LINE!

	silent_data.close()

	COL_NAME_LIST=COLUMN_NAME_LINE.split()

	if(SEQUENCE_LINE[0:9]!='SEQUENCE:'):  		 error_exit_with_message("SEQUENCE_LINE[0:9]!='SEQUENCE:' for SEQUENCE_LINE (%s)" %(SEQUENCE_LINE) )
	if(COLUMN_NAME_LINE[0:6]!='SCORE:'):  		 error_exit_with_message("COLUMN_NAME_LINE[0:6]!='SCORE:'" )
	if(COLUMN_NAME_LINE.find('description') == -1 ):  error_exit_with_message("COLUMN_NAME_LINE.find('description') == -1" )
	if(COL_NAME_LIST[0]!='SCORE:'): error_exit_with_message("COL_NAME_LIST[0]!='SCORE:" )
	if(COL_NAME_LIST[-1]!='description'): error_exit_with_message("COL_NAME_LIST[-1]!='description'" )

	##Here we will assume that


	#firstlines = popen_and_readlines('head -n 3 '+ silent_file, Is_master=False, tag=("head_" + silent_file.replace('/','_') ) )

	try:
		score_col_index=COL_NAME_LIST.index(scorecol_name)
	except:
		error_exit_with_message("Cannot find score_col_index with corresponding score_colname= %s " %(scorecol_name) )

	try:
		tag_col_index=COL_NAME_LIST.index('description')
	except:
		error_exit_with_message("Cannot find description column index! ")

	print "%s is located at column_num: %d" %(scorecol_name, score_col_index+1),
	print "tag is located at column_num: %d" %(tag_col_index+1)
	
	################################################################################################################################################################
	  
	# Get the SCORE line    
	SCORE_silentfile=dirname(silent_file) + '/SCORE_' + basename(silent_file)

	if(SCORE_silentfile[0]=='/'): SCORE_silentfile=SCORE_silentfile[1:] #Oct 3, 2010

	command = 'grep "SCORE: " %s > %s' %(silent_file, SCORE_silentfile)
	if(verbose): print command
	sys.stdout.flush()
	submit_subprocess(command)

	############################Get the extremum_score########################################
	if(verbose): print "Get the extremum_score"	
	sys.stdout.flush()

	extremum_score=0

	if(mode_score_gap == True):

		try:
			infile = open( SCORE_silentfile, 'r' )
		except:
			error_exit_with_message("cannot open %s " % SCORE_silentfile)

		first_line=True
		for line in infile:
		
			if(Is_column_name_line(line)==True): continue ##Take this off it is significantly slow down code
		
			cols = string.split( line )
		
			try:
				score = float( cols[ score_col_index ] )
			except:
				print "problem with reading score line"
				print "SCORE_silentfile: %s" %(SCORE_silentfile)
				print "score_col_index= ", score_col_index 
				print "cols: ", cols
				error_exit_with_message("problem with reading score line")
		
			if(first_line==True): 
				extremum_score=score
				first_line=False
		
			if(REVERSE == True): #extremum is a minimum
				if(extremum_score > score): extremum_score=score 	
			else:                #extremum is a maximum
				if(extremum_score < score): extremum_score=score
			
		infile.close()
		
	if(verbose): print "exteremum_score= %f" %(extremum_score)	
	############################Create the score_and_tag_list#################################
	if(verbose): print "Create the score_and_tag_list"	
	sys.stdout.flush()


	try:
		infile = open( SCORE_silentfile,'r')
	except:
		error_exit_with_message("cannot open %s " %(SCORE_silentfile))
	
	score_and_tag_list=[]

	for line in infile:
	
		if(Is_column_name_line(line)==True): continue ##Take this off it is significantly slow down code

		cols = string.split( line )

		tag=cols[tag_col_index]

		try:
			score = float( cols[ score_col_index ] )
		except:
			print "problem with reading score line"
			print "SCORE_silentfile: %s" %SCORE_silentfile
			print "score_col_index= ", score_col_index 
			print "cols: ", cols
			error_exit_with_message("problem with reading score line")
	
		###Check that within E_gap here...
		if(mode_score_gap == True):
			if(REVERSE==True): #extremum is a minimum
				if(extremum_score+max_score_gap < score): continue
			else:
				if(extremum_score-max_score_gap > score): continue


		if(mode_abs_score_cut == True):
			if(REVERSE==True): #extremum is a minimum
				if(abs_score_cut < score): continue
			else:
				if(abs_score_cut > score): continue
			
		
		score_and_tag_list.append( ( score, tag ) ) 

	infile.close()

	##############Last update on Jan 14, 2012###########################################
	if(exists(SCORE_silentfile)==False): error_exit_with_message("SCORE_silentfile (%s) doesn't exist!" %(SCORE_silentfile) )
	if(remove_SCORE_file):	
		print "removing SCORE_silentfile (%s) ......" %(SCORE_silentfile)
		submit_subprocess("rm %s" %(SCORE_silentfile))  
	####################################################################################

	###So by default it sort by the 0th element which is what we want
	###If reverse==false, then it sort from low to high which is what we want typically..if reverse==true then it sort from high to low
#	score_and_tag_list.sort(reverse=(not REVERSE))

	###sort() will sort the list from low to high which is what we want typically
	score_and_tag_list.sort()

	#if REVERSE=false, then reverse the list since we now want if from high to low
	if(REVERSE==False): score_and_tag_list.reverse() 
 		 

	score_and_tag_list=score_and_tag_list[:int(max_n_struct)]	###Only keep the best NSTRUCT datapoints ##NSTRUCT will be infinite in mode_score_gap. 

	actual_num_structures=len(score_and_tag_list)
	actual_score_gap=score_and_tag_list[0][0]-score_and_tag_list[len(score_and_tag_list)-1][0]

	##Have a if statement for consistency check here

	if(verbose): print "actual_num_structures= %d actual_score_gap=%f (inside SWA_filter_silent_file) "	%(actual_num_structures, actual_score_gap)
	sys.stdout.flush()

	#####################Create the list of tags that will be kept#######################################
	tags = []

	for score_and_tag in score_and_tag_list:
		tag = score_and_tag[1] ##2nd element
		tags.append(tag)

	tags_set=Set(tags) #convert to a set Sept 28, 2010...membership test is much faster this way!
	
	####################Create the outfile##################################
	print "start writing to filter_outfile %s " %(outfile_name)
	start_writing_time=time.time()

	try:
		outfile = open( outfile_name , 'w')  # this removes the old outfile_name if one exist, for example when refiltering the silent_file
	except:
		error_exit_with_message("cannot open outfile (%s) " %(outfile))

	outfile.write(SEQUENCE_LINE) #Sequence line
	outfile.write(COLUMN_NAME_LINE) #Column name line
	if(THIRD_LINE[:6] == 'REMARK' ): outfile.write(THIRD_LINE) #REMARK line

	try:
		infile = open( silent_file , 'r')
	except:
		error_exit_with_message("cannot open infile (%s) " %(infile))

	output_count=0
	next_output_time=0.0
		
	writeout = False
	for line in infile:
		cols = string.split(line)
		if(len(cols)>1 and cols[0]=='SCORE:'):

			if( cols[tag_col_index] in tags_set ):
				writeout = True
				output_count+=1
			else:
				writeout = False
	
		#########################		
		if(writeout): 
			outfile.write(line)

			time_taken=time.time()-start_writing_time
			if( time_taken>next_output_time ):
				next_output_time+=2.0 #seconds
				print "%d struct output/written in %f seconds" %(output_count, time_taken) 	
				sys.stdout.flush()
				sys.stderr.flush()	
		#########################

	outfile.close()
	infile.close()

	time_taken=time.time()-start_writing_time
	print "FINAL: %d struct output/written in %f seconds" %(output_count, time_taken) 			
	print_title_text("Exit " + title)
	sys.stdout.flush()
	sys.stderr.flush()

	return (actual_num_structures, actual_score_gap)

