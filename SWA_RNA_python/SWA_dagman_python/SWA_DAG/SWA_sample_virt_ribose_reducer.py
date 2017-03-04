#!/usr/bin/env python

from sys import argv
import sys
import traceback
from glob import glob
import string
from os.path import basename, dirname, exists, expanduser
from time import sleep
import time
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################


#############################################################################

#SWA_sample_virt_ribose_reducer.py -data_foldername VIRTUAL_RIBOSE_SAMPLER/REGION_0_2_FOR_REGION_3_2_START_FROM_REGION_0_2_AND_4_0/ -reducer_outfile reducer_sample_virt_ribose_region_0_2.out -START_silent_file region_0_2_sample.cluster.out -test_mode True -keep_old_score True -keep_old_tag True

python_command=get_python_command(argv)
print_title_text("Enter " + python_command)

data_foldername = parse_options( argv, "data_foldername", "" )

if(data_foldername == ""): error_exit_with_message("data_foldername == \"\"!")
if(exists(data_foldername)==False): error_exit_with_message("data_foldername (%s) doesn't exist!" %(data_foldername) )


reducer_outfile = parse_options( argv, "reducer_outfile", "" )

if(reducer_outfile == ""): error_exit_with_message("reducer_outfile == \"\"!")
if(dirname(reducer_outfile) != ""):  error_exit_with_message("dirname(reducer_outfile) != \"\"")

test_mode=parse_options( argv, "test_mode", "False" )
delete_files=True

if(test_mode==False): reducer_outfile=data_foldername + '/' + reducer_outfile

if(test_mode): delete_files=False

print "test_mode=%s | delete_files=%s" %(test_mode, delete_files)

START_silent_file = parse_options( argv, "START_silent_file", "" )

if( START_silent_file == "" ): error_exit_with_message("START_silent_file == \"\"!")
if(exists(START_silent_file)==False): error_exit_with_message("START_silent_file (%s) doesn't exist!" %(START_silent_file) )

keep_old_score= parse_options( argv, "keep_old_score", "False" )
keep_old_tag= parse_options( argv, "keep_old_tag", "False" )

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )

########################################################################################################################

deleting_files_signal_file="%s/SWA_sample_virt_ribose_reducer_deleting_files_signal_file.txt" %(data_foldername)

if(exists(deleting_files_signal_file)):	error_exit_with_message("%s exist!" %(deleting_files_signal_file))

#if(exists(reducer_outfile)==True): error_exit_with_message("reducer_outfile (%s) already exist!" %(reducer_outfile) ) #COMMENTED OUT ON FEB 09, 2012

if( exists(reducer_outfile) ): #FEB 09, 2012
	print "Warning reducer_outfile (%s) already exist....removing!" %(reducer_outfile)
	submit_subprocess("rm %s" %(reducer_outfile) )

########################################################################################################################

### check for valid non_empty_silent_file
if ( not Is_valid_non_empty_silent_file( START_silent_file ) ):
	
	### EARLY EXIT if valid_empty_silent_file ... 
	if ( Is_valid_empty_silent_file( START_silent_file ) ):
		print "START_silent_file (%s) is a valid_empty_silent_file!" % (START_silent_file)
		with open( reducer_outfile , 'w' ) as REDUCER_OUTFILE:
			REDUCER_OUTFILE.write( "empty cluster silent_file since all input_silent_file are empty." )
		print_title_text("Exit " + python_command)
		sys.exit(0)

	else:
		error_exit_with_message("START_silent_file (%s) is not a valid_non_empty_silent_file!" %( START_silent_file ) )
	
########################################################################################################################

num_column_name_line=0
#count=0 ##Commented out on Feb 08, 2012

data = safe_open(START_silent_file, mode='r', Is_master=False)
SEQUENCE_LINE = data.readline()
COLUMN_NAME_LINE = data.readline()
THIRD_LINE = data.readline() #Possible REMARK LINE!
third_line_is_a_remark = ("REMARK" in THIRD_LINE)
data.close()
print "SEQUENCE_LINE=", SEQUENCE_LINE.replace('\n', '')
print "COLUMN_NAME_LINE=", COLUMN_NAME_LINE.replace('\n', '')
print "THIRD_LINE=", THIRD_LINE.replace('\n', '')
print "third_line_is_a_remark=", str(third_line_is_a_remark)

COL_NAME_LIST = COLUMN_NAME_LINE.split()

try:
	score_col_index = COL_NAME_LIST.index('score')
except:
	print "COL_NAME_LIST=", COL_NAME_LIST
	error_exit_with_message("Cannot find score column index!" )

try:
	tag_col_index=COL_NAME_LIST.index('description')
except:
	print "COL_NAME_LIST=", COL_NAME_LIST
	error_exit_with_message("Cannot find description column index!")

START_silent_data_list = []

offset = 0

data = safe_open(START_silent_file, mode='r', Is_master=False)
dummy_sequence_line = data.readline()
dummy_column_name_line = data.readline()
if third_line_is_a_remark:	
	dummy_third_line = data.readline() # only read third line if it is remark; otherwise, it could be the SCORE line.

while(True):

	offset=data.tell()

	line=data.readline()

	if(line==''): break #End of file!

	if(len(line) <= 1): error_exit_with_message("len(line) <= 1") #check for line with only '\n'

	if( (line.find("SCORE") != -1) != (line[0:6]=="SCORE:") ): error_exit_with_message("(line.find(\"SCORE\") != -1) != (line[0:6]==\"SCORE:\")")

	if(line[0:6]=="SCORE:"): #SCORE line!

		if(line.find('description') != -1): error_exit_with_message("extra column_name line (%s)" %(line) )

		score_line=line

		#count+=1 ##Commented out on Feb 08, 2012
		#if(test_mode and count>5): break ##Commented out on Feb 08, 2012

		tag=score_line.split()[tag_col_index]

		score_string=score_line.split()[score_col_index]
		score=float(score_string)

		queue_ID=int(tag.split("_")[1])
		if(queue_ID!=len(START_silent_data_list) ): error_exit_with_message("queue_ID(%d)!=len(START_silent_data_list)(%d)" %(queue_ID, len(START_silent_data_list) ) )

		for n in range(len(START_silent_data_list) ):
			if(START_silent_data_list[n]["tag"]==tag): error_exit_with_message("silent_data with tag (%s) already exist in the START_silent_data_list!" %(tag))

		silent_data={}
		silent_data["tag"]=tag
		silent_data["queue_ID"]=queue_ID
		silent_data["score_string"]=score_string
		silent_data["score"]=score
		silent_data["file_name"]=START_silent_file
		silent_data["offset"]=offset

		START_silent_data_list.append(silent_data)


	offset=data.tell()

data.close()

if(len(START_silent_data_list)==0): error_exit_with_message("len(START_silent_data_list)==0")

globstring = data_foldername+'/S_*/' + 'sample_virt_ribose_region_*.out'
print "globstring= ", globstring
sys.stdout.flush()

globfiles = glob( globstring )
if(len( globfiles ) == 0):
	sleep( 10 )
	globfiles = glob( globstring )
globfiles.sort()

globfile_info_list=[]

for n in range(len(globfiles)):

	globfile_info={}

	globfile_info["file_name"]=globfiles[n]
	globfile_info["tag"]=globfile_info["file_name"].split('/')[-2]
	globfile_info["queue_ID"]=int(globfile_info["tag"].split("_")[-1])

	globfile_info_list.append(globfile_info)

globfile_info_list=sort_dictionary_list(globfile_info_list, "queue_ID")

if(len(globfile_info_list)!=len(START_silent_data_list)):
	print "len(globfile_info_list)=%d" %(len(globfile_info_list))
	print "len(START_silent_data_list)=%d" %(len(START_silent_data_list))
	error_exit_with_message("len(globfile_info_list)!=len(START_silent_data_list)")

NEW_silent_data_list=[]

num_fail_ribose_sampling=0

for n in range(len(START_silent_data_list)):

	START_silent_data=START_silent_data_list[n]
	globfile_info    =globfile_info_list[n]

	###consistency checks
	if(n!=START_silent_data["queue_ID"]): error_exit_with_message("n!=START_silent_data[\"queue_ID\"] for n=%d " %(n) )

	if(n!=globfile_info["queue_ID"]): error_exit_with_message("n!=globfile_info_list[\"queue_ID\"] for n=%d" %(n) )

	if(globfile_info["tag"]!=START_silent_data["tag"]): error_exit_with_message( "globfile_info[\"tag\"]!=START_silent_data[\"tag\"]" )

	data=safe_open(globfile_info["file_name"], mode='r', Is_master=False)

	first_line=data.readline();

	data.close()

	num_options_taken=0

	if(first_line=="no_virtual_ribose (FB_CC_JP_list.size()==0).\n") or \
	  (first_line=="no_virtual_sugar ( sugar_modeling_list.size() == 0 ).\n" ): #No virtual_ribose
	
		NEW_silent_data_list.append(START_silent_data)

	elif(first_line=="num_virtual_ribose != 0 but for one of the sampled virtual_ribose, curr_FB_JP.PDL.size()==0.\n" ) or \
		(first_line=="num_virtual_sugar != 0 but for one of the sampled virtual_sugar, curr_sugar_list.size()==0.\n" ) or \
		(first_line=="num_virtual_sugar != 0 but for one of the sampled virtual_sugar,curr_modeling.pose_list.size() == 0.\n" ): #Virtual_ribose exist, but no valiable ribose rotamer!
	
		num_fail_ribose_sampling+=1 # Filler to prevent error!

	else:

		assert_is_valid_non_empty_silent_file(globfile_info["file_name"])

		data = safe_open(globfile_info["file_name"], mode='r', Is_master=False)

		first_line = data.readline()
		second_line = data.readline()
		third_line = data.readline()

		error1 = (SEQUENCE_LINE != first_line)
		error2 = (not Is_equivalent_list(COL_NAME_LIST, second_line.split()))
		error3 = ((THIRD_LINE[:6] == 'REMARK') and (THIRD_LINE != third_line))

		if error1 or error2 or error3:	
			
			print "\nfile_name=%s" % (globfile_info["file_name"])

			if error1:	#(SEQUENCE_LINE!=first_line):
				print "SEQUENCE_LINE =", SEQUENCE_LINE.replace('\n','')
				print "first_line    =", first_line.replace('\n','')
				### If full_model_info is defined, the full_model sequence will be written to the silent file, causing some reducing steps to error out.
				#error_exit_with_message("SEQUENCE_LINE !=first_line")

			if error2:	#( Is_equivalent_list(COL_NAME_LIST, second_line.split())==False ) :
				print "COL_NAME_LIST      =", COL_NAME_LIST
				print "second_line.split()=", second_line.split()
				error_exit_with_message("Is_equivalent_list(COL_NAME_LIST, second_line.split()==False)" )

			if error3:	#( (THIRD_LINE[:6] == 'REMARK') and (THIRD_LINE!=third_line) ):
				print "THIRD_LINE=", THIRD_LINE.replace('\n','')
				print "third_line=", third_line.replace('\n','')
				#error_exit_with_message("(THIRD_LINE[:6] == 'REMARK') and (THIRD_LINE!=third_line))" )

		offset=0

		while(True):

			offset=data.tell()

			line=data.readline()

			if(line==''): break #End of file!

			if(len(line) <= 1): error_exit_with_message("len(line) <= 1") #check for line with only '\n'

			if( (line.find("SCORE") != -1) != (line[0:6]=="SCORE:") ): error_exit_with_message("(line.find(\"SCORE\") != -1) != (line[0:6]==\"SCORE:\")")

			if(line[0:6]=="SCORE:"): #SCORE line!

				if(line.find('description') != -1): error_exit_with_message("extra column_name line (%s)" %(line) )


				score_line=line

				#count+=1 ##Commented out on Feb 08, 2012
				#if(count==10): break ##Commented out on Feb 08, 2012
				#THIS IS AN ERROR..BUT HOPEFULLY DOES NOT SIGNIFICANTLY EFFECT RESULTS SINCE
				#1. If silent_file have more than 10 structures, the count is greater than 10 before reaching this break statement(see above)
				#2. AT IN WORSE CASE if reach this break statement, ONLY 1 STRUCTURE (possibly when multiple sibling ribose conformations) is IGNORED.

				new_silent_data={}

				new_silent_data["file_name"]=globfile_info["file_name"]

				new_silent_data["tag"]=score_line.split()[tag_col_index]

				if(keep_old_score):
					new_silent_data["score_string"]=START_silent_data["score_string"]
					new_silent_data["score"]=START_silent_data["score"]
				else:
					new_silent_data["score_string"]=(score_line.split()[score_col_index])
					new_silent_data["score"]=float(new_silent_data["score_string"])

				new_silent_data["queue_ID"]=int(new_silent_data["tag"].split("_")[1])

				new_silent_data["offset"]=offset

				parent_tag="S_%d" %(new_silent_data["queue_ID"])

				if(parent_tag!=START_silent_data["tag"]): error_exit_with_message("parent_tag!=START_silent_data[\"tag\"]") #Consistency check.

				NEW_silent_data_list.append(new_silent_data)

		data.close()

print "num_fail_ribose_sampling=%d" %(num_fail_ribose_sampling)

if(keep_old_score==False): #Resort the silent_file_data!

	NEW_silent_data_list=sort_dictionary_list(NEW_silent_data_list, key="score")


for n in range( len(NEW_silent_data_list) ):

	if(keep_old_tag):
		NEW_silent_data_list[n]["new_tag"]=NEW_silent_data_list[n]["tag"]
	else:
		NEW_silent_data_list[n]["new_tag"]="S_%d" %(n)

	NEW_silent_data_list[n]["new_queue_ID"]= n

###Important to assert that the new_tags are unique###
new_tag_set=Set([]) #empty set!

for n in range( len(NEW_silent_data_list) ):

	new_tag=NEW_silent_data_list[n]["new_tag"]

	if(new_tag in new_tag_set): error_exit_with_message("new_tag (%s) already exist in new_tag_set!")

	new_tag_set.add(new_tag)

new_tag_set.clear()

if(test_mode):
	for n in range(len(NEW_silent_data_list)):
		print NEW_silent_data_list[n]


print "--------------------start writing to reducer_outfile %s--------------------" %(reducer_outfile)
start_writing_time=time.time()

REDUCER_OUTFILE= open( reducer_outfile , 'w')  # this removes the old outfile_name if one exist

REDUCER_OUTFILE.write(SEQUENCE_LINE) #Sequence line
REDUCER_OUTFILE.write(COLUMN_NAME_LINE) #Column name line ####DOES it matter that the column width of the score_line for different silent_file doesn't match?
if( (THIRD_LINE[:6] == 'REMARK') ): REDUCER_OUTFILE.write(THIRD_LINE)

output_count=0
next_output_time=0.0

for n in range(len(NEW_silent_data_list)):

	output_count+=1

	silent_data=NEW_silent_data_list[n]

	infile=safe_open(silent_data["file_name"], mode='r', Is_master=False)

	writeout = False

	infile.seek(silent_data["offset"])

	score_line=infile.readline()

	cols = score_line.split()

	if( len(cols) != len(COL_NAME_LIST)): error_exit_with_message("len(cols) != len(COL_NAME_LIST) ")

	if( len(cols) <= tag_col_index ): error_exit_with_message("len(cols)<=tag_col_index")

	if( cols[tag_col_index] != silent_data["tag"] ): error_exit_with_message("cols[tag_col_index] != silent_data[\"tag\"]")

	new_cols = cols

	if(keep_old_score):
		new_cols[score_col_index] = silent_data["score_string"] #we have replaced score_string with "old" score_string above
	else:
		#consistency_check
		if(new_cols[score_col_index]!= silent_data["score_string"]): error_exit_with_message("new_cols[score_col_index]!= silent_data[\"score_string\"]")

	new_cols[tag_col_index] = silent_data["new_tag"]

	REDUCER_OUTFILE.write( list_to_string(new_cols, first_separator=False) +'\n' ) #Ok this will change the blank spacing between the column. Does this matter?

	#REDUCER_OUTFILE.write( score_line[:score_line.rindex(' ')+1 ] + silent_data["new_tag"]+'\n' )

	REDUCER_OUTFILE.write("REMARK PARENT_TAG %s\n" %(silent_data["tag"] ) )
	REDUCER_OUTFILE.write("REMARK SOURCE %s\n" %(silent_data["file_name"]) )

	for line in infile:

		# Found the next score_line! break!
		if(line.split()[0]=='SCORE:'): break

		if(line[0:17]=="REMARK PARENT_TAG"): continue
		if(line[0:13]=="REMARK SOURCE"): continue

		if(line[0:7]=="REMARK "):	 #other remarks?
			REDUCER_OUTFILE.write( line )
		else:
			REDUCER_OUTFILE.write( line[:line.rindex(' ')+1 ] + silent_data["new_tag"]+'\n' )


	time_taken=time.time()-start_writing_time

	if( time_taken>next_output_time ):
		next_output_time+=2.0 #seconds
		print "%d struct output/written in %f seconds" %(output_count, time_taken)
		sys.stdout.flush()
		sys.stderr.flush()

	infile.close()

time_taken=time.time()-start_writing_time
print "FINAL: %d struct output/written in %f seconds" %(output_count, time_taken)
sys.stdout.flush()
sys.stderr.flush()


REDUCER_OUTFILE.close()


if(delete_files):

	create_generic_done_signal_file(deleting_files_signal_file)		#The mark the point where the scipt is now not rerunnable!

	for filename in globfiles:
		foldername=dirname(filename)
		if(foldername==""): error_exit_with_message("foldername==\"\"")
		delete_command="rm -r %s " %(foldername)
		print delete_command
		sys.stdout.flush()
		sys.stderr.flush()
		submit_subprocess( delete_command )

print_title_text("Exit " + python_command)

