#!/usr/bin/env python

import string
from sys import argv
from os.path import exists,basename

######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.utility.DAGMAN_util import update_CONDOR_file_with_actual_job_queue_num

######################################################################


python_command = get_python_command(argv)
print_title_text("Enter: " + python_command)


filterer_outfile = argv[1]  
condor_submit_file = argv[2] 

#CHECK if the input silent_file is empty###############
#firstline = popen_and_readlines('head -n 1 '+ filterer_outfile, Is_master=False, tag=("head_1_" + filterer_outfile.replace('/','_') + '_' + condor_submit_file.replace('/','_') ) )[0]
#found_early_exit_signal=False

#if( firstline.count("Empty filterer_outfile. No struct_pair passed screen.")!=0):
#	found_early_exit_signal=True

#if( firstline.count("empty cluster silent_file since at least one of the two input_silent_file is empty.")!=0):
#	found_early_exit_signal=True

#if( firstline.count("empty cluster silent_file since at least of the two input_silent_file is empty.")!=0): #this one is now obsolete
#	found_early_exit_signal=True

#if( found_early_exit_signal ):
if ( Is_valid_empty_silent_file( filterer_outfile ) ):

	print "EARLY EXIT: filterer_outfile %s contain no tags pair " %(filterer_outfile)

	update_CONDOR_file_with_actual_job_queue_num(condor_submit_file, 0)
	print_title_text("Exit: " + python_command) #Early exit
	sys.exit(0)


#################################################################################################################


try:
	infile = open( filterer_outfile ,'r')
except:
	error_exit_with_message("cannot open %s " % filterer_outfile)
	
job_num=-1

for line in infile:
	
#	print line

	line_list=line.split()

	append_side_tag=line_list[0]
	prepend_side_tag=line_list[1]

	job_num+=1


	print 'append_side_tag= %s, prepend_side_tag= %s' %(append_side_tag, prepend_side_tag)

infile.close()



##################################################################################################################
N_JOBS = job_num+1

update_CONDOR_file_with_actual_job_queue_num(condor_submit_file, N_JOBS)

print_title_text("Exit: " + python_command)
