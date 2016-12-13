#!/usr/bin/env python

from sys import argv,exit,stderr
from os import system
from os.path import basename
import string
from SWA_util import *

#print "max_score_gap=", max_score_gap
def filter_silent_file(silent_file, outfile_name, scorecol_name="score", REVERSE=True, max_n_struct=0, max_score_gap=0.0):

	title='filter_silent_file()'

	print_title_text("Enter " + title)

	print "silent_file=%s " %silent_file
	print "outfile_name=%s " %outfile_name
	print "REVERSE= ", REVERSE

	if(silent_file == ""):
		error_exit_with_message("need to input silent_file!")
	if(outfile_name == ""):
		error_exit_with_message("need to outfile_name")


	mode_max_struct=False
	mode_score_gap=False

	print "scorecol_name=%s" %(scorecol_name)

	if(max_n_struct !=0): mode_max_struct=True
	if(max_score_gap > 0.01): mode_score_gap=True

	if( mode_max_struct==False and mode_score_gap==False):
		error_exit_with_message("need to specify n_struct or score_gap")
	elif(mode_max_struct==True and mode_score_gap==True):
		error_exit_with_message("Cannot specify both n_struct and score_gap ")
	elif(mode_max_struct):
		print "mode_max_struct, max_n_struct=%d, max_score_gap=%8.3f" %(max_n_struct, max_score_gap)
	elif(mode_score_gap):
		max_n_struct=99999999999999999999
		if(max_score_gap<0.0): error_exit_with_message("max_score_gap= %8.3f need to be a positive real number" %(max_score_gap))
		print "mode_score_gap, max_n_struct=%d, max_score_gap=%8.3f" %(max_n_struct, max_score_gap)
	else:
		error_exit_with_message("Should not reach this point of the code! ")
	#######################################################################################################################################
	##Here we will assume that
	##First line: SEQUENCE
	##Second line: Columns name
	##Third line: possible a REMARK LINE or else start of normal data

	firstlines = popen_and_readlines('head -n 3 '+ silent_file)

	col_name_line=firstlines[1]
	col_name_list=string.split(col_name_line)

	try:
		score_col_index=col_name_list.index(scorecol_name)
	except:
		error_exit_with_message("Cannot find score_col_index with corresponding score_colname= %s " %(scorecol_name) )

	try:
		tag_col_index=col_name_list.index('description')
	except:
		error_exit_with_message("Cannot find description column index! ")

	print "%s is located at column_num: %d" %(scorecol_name, score_col_index+1),
	print "tag is located at column_num: %d" %(tag_col_index+1)

	################################################################################################################################################################



	# Get the SCORE line
	SCORE_silentfile=dirname(silent_file) + '/SCORE_' + basename(silent_file)

#	command = 'grep SCORE %s > %s' %(silent_file, SCORE_silentfile)  buggy!
	command = 'grep "SCORE: " %s > %s' %(silent_file, SCORE_silentfile)

	submit_subprocess(command)

	############################Get the extremum_score########################################

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
				print "SCORE_silentfile: %s" %SCORE_silentfile
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

	############################Create the score_and_tag_list#################################
	try:
		infile = open( SCORE_silentfile,'r')
	except:
		error_exit_with_message("cannot open %s " % SCORE_silentfile)

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

		score_and_tag_list.append( ( score, tag ) )

	infile.close()
#	submit_subprocess("rm %s" %SCORE_silentfile) #Mod this: Don't remove SCORE_silentfile. Keep it for consistency, doesn't take up too much space anyways.


	###So by default it sort by the 0th element which is what we want
	###If reverse==false, then it sort from low to high which is what we want typically..if reverse==true then it sort from high to low
#	score_and_tag_list.sort(reverse=(not REVERSE))

	###sort() will sort the list from low to high which is what we want typically
	score_and_tag_list.sort()

	if(not REVERSE): #if REVERSE=false, then reverse the list since we now want if from high to low
 		score_and_tag_list.reverse()


	score_and_tag_list=score_and_tag_list[:int(max_n_struct)]	###Only keep the best NSTRUCT datapoints ##NSTRUCT will be infinite in mode_score_gap.

	actual_num_structures=len(score_and_tag_list)
	actual_score_gap=score_and_tag_list[0][0]-score_and_tag_list[len(score_and_tag_list)-1][0]

	##Have a if statement for consistency check here

	#####################Create the list of tags that will be kept#######################################
	tags = []

	templist_name = 'temp.%s.list' %(basename(silent_file))  ####A list of tags

	fid = open(templist_name,'w')

	for score_and_tag in score_and_tag_list:
		tag = score_and_tag[1] ##2nd element
		tags.append(tag)

		####Debug output
		fid.write('%s	%f'"\n" %(tag, score_and_tag[0]))


	fid.close()

	command = 'rm '+templist_name
	print(command)
	submit_subprocess(command) #Mod this for testing for now

	####################Create the outfile##################################
	outfile = open( outfile_name , 'w')

	outfile.write(firstlines[0]) #Sequence line
	outfile.write(firstlines[1]) #Column name line

	if (firstlines[2][:6] == 'REMARK' ):
		outfile.write(firstlines[2]) #REMARK line

	infile = open( silent_file , 'r')


	writeout = False
	for line in infile:
		cols = string.split(line)
		if (len(cols)>1 and cols[0]=='SCORE:'):
			if tags.count(cols[tag_col_index]) > 0:
				writeout = True
			else:
				writeout = False
		if writeout:
			outfile.write(line)

	outfile.close()
	infile.close()

	print_title_text("Exit " + title)

	return (actual_num_structures, actual_score_gap)

###################################################Run from command line########################################################################
'''
if(len(argv)>1):

	python_command=get_python_command(argv)
	print_title_text("Enter " + get_python_command(argv))

	Parse_silent_file = parse_options( argv, "silent_file", "" )
	if(silent_file == ""):
		error_exit_with_message("need to input silent_file!")

	Parse_outfile_name = parse_options( argv, "outfile_name", "" )
	if(outfile_name == ""):
		error_exit_with_message("need to outfile_name")

	#Python Boolean. If True then sort from low to high. If False then from high to low
	Parse_REVERSE = parse_options( argv, "reverse", "True" )

	##Basically scorecol is used in a general sense here. scorecol Can be any column: rmsd, score and etc.
	Parse_scorecol_name = parse_options( argv, "scorecol_name", "score" )

	Parse_max_n_struct=parse_options( argv, "n_struct", 0 ) ##Num struct kept
	MAR 25, WARNING max_score_gap is now a REAL NUMBER
	Parse_max_score_gap=parse_options( argv, "score_gap", 0)
	MAR 25, WARNING max_score_gap is now a REAL NUMBER
	filter_silent_file(Parse_silent_file, Parse_outfile_name, scorecol_name=Parse_scorecol_name, REVERSE=Parse_REVERSE, max_n_struct=Parse_max_n_struct, max_score_gap=Parse_max_score_gap)
	print_title_text("Exit " + python_command)
'''
