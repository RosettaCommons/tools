#!/usr/bin/env python

import string
from sys import argv
from os.path import exists,basename
from random import random

######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.utility.DAGMAN_util import update_CONDOR_file_with_actual_job_queue_num

######################################################################

###Example: /home/sripakpa/SWA_dagman_python/SWA_sampling_pre_process.py REGION_5_1 REGION_5_0 CONDOR/REGION_5_1_START_FROM_REGION_5_0.condor

# random delay -- needed for condor_dagman which sometimes launches pre_process scripts at exactly the same time which try to read the same files!
print "Random delay of about 2 seconds..."
sleep( 2 * random() )

python_command = get_python_command(argv)
print_title_text("Enter: " + python_command)

input_silent_file  = argv[1]
condor_submit_file = argv[2]

#check that there is the string input_silent_file is a single filename
if(len(input_silent_file.split())!=1):
	error_exit_with_message("input_silent_file.split()!=1, input_silent_file= %s" %(input_silent_file) )

#input_silent_file='%s_sample.cluster.out' %(dir_prev.lower())  ###AGAIN this is hard coded

#CHECK if the input silent_file is empty###############
#first_line = popen_and_readlines('head -n 1 '+ input_silent_file, Is_master=False, tag=("head_1_" + input_silent_file.replace('/','_') + '_' + condor_submit_file.replace('/','_') ) )[0]
#print "first_line= ", first_line
#if ( "empty cluster silent_file since all input_silent_file are empty." in first_line ):
if ( Is_valid_empty_silent_file( input_silent_file ) ):

	print "EARLY EXIT: input_silent_file %s contain no silent structs " % (input_silent_file)

	update_CONDOR_file_with_actual_job_queue_num(condor_submit_file, 0)
	print_title_text("Exit: " + python_command) #Early exit
	sys.exit(0)


#################################################################################################################
#run by master node, but through a submit_process()
firstlines = popen_and_readlines('head -n 3 '+ input_silent_file, Is_master=False, tag=("head_3_" + input_silent_file.replace('/','_') + '_' + condor_submit_file.replace('/','_') ) )

col_name_line=firstlines[1] ###AGAIN this assumes that the second line of the silent file is the column name line
col_name_list=string.split(col_name_line)

try:
	tag_col_index=col_name_list.index('description')
except:
	error_exit_with_message("Cannot find description column index! ")

print "tag is located at column_num: %d" %(tag_col_index+1)

SCORE_silentfile='SCORE_' + basename(input_silent_file) + '_' +  condor_submit_file.replace('/','_')

command = 'grep "SCORE: " %s > %s' %(input_silent_file, SCORE_silentfile) #Sept 14.. the old version is buggy since the word SCORE can appear in the binary struct code

#command = 'grep SCORE %s > %s' %(input_silent_file, SCORE_silentfile)

submit_subprocess_allow_retry(command)

try:
	infile = open( SCORE_silentfile,'r')
except:
	error_exit_with_message("cannot open %s " % SCORE_silentfile)

###Note: Elsewhere, the python code had been hard coded to assume that the tags from input_silentfile is of the form:
### S_0, S_1, S_2, S_3,...... Personally I think we should relax this assumption but  to do thiswill have to make
### many changes elsewhere in the code. The best I could do for now is to ensure here that the input_silentfile indeed
### does have this specific form and raise an error otherwise.

pose_num=-1

for line in infile:

	if(Is_column_name_line(line)==True): continue

	pose_num+=1
	assumed_tag = 'S_%d' %pose_num
#	assumed_tag = 'S_%06d' %pose_num

#S_000000

	cols = string.split( line )
#	print cols

	actual_tag=cols[tag_col_index]

	if(assumed_tag!=actual_tag):
		error_exit_with_message("assumed_tag (%s)!=actual_tag (%s)" %(assumed_tag, actual_tag) )



infile.close()
submit_subprocess_allow_retry("rm %s" %SCORE_silentfile)


##################################################################################################################
N_JOBS = pose_num+1

update_CONDOR_file_with_actual_job_queue_num(condor_submit_file, N_JOBS)

print_title_text("Exit: " + python_command)

sys.exit( 0 )
