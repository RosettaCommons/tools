#!/usr/bin/env python


from os import system,popen
from os.path import exists,dirname,basename,expanduser,abspath
from sys import exit, argv
import string
from time import sleep
import glob
import time
import fileinput
import os
######################################################################

from SWA_dagman_python.misc.SWA_cat_outfiles import concatenate_outfiles
######################################################################


from SWA_dagman_python.parser.SWA_parse_options import parse_options
from SWA_dagman_python.utility.SWA_util import *
######################################################################

#Examples: 

#DAG_rebuild_bulge_reducer.py -outfolder FINAL_REBUILD_BULGE/STANDARD/ -reducer_outfile rebuild_bulge_region_FINAL.out -condor_submit_file CONDOR/REBUILD_BULGE/FINAL_REBUILD_BULGE.condor  

#DAG_rebuild_bulge_reducer.py -outfolder TEST/ -reducer_outfile rebuild_bulge_region_FINAL.out -condor_submit_file MOCK_FINAL_REBUILD_BULGE.condor

#run by master node, but through a submit_process()

python_command = get_python_command(argv)
print_title_text("Enter: " + python_command)

data_foldername= parse_options( argv, "outfolder", "" ) #FINAL_REBUILD_BULGE/STANDARD/

if(data_foldername==""): error_exit_with_message("User need to pass in data_foldername ")

if(exists( data_foldername )==False): error_exit_with_message("data_foldername (%s) doesn't exist!" %(abspath(data_foldername) ))

condor_submit_file=parse_options( argv, "condor_submit_file", "" ) #FEB 08, 2012
if(condor_submit_file == ""): error_exit_with_message("condor_submit_file == \"\" ") 

reducer_outfile= parse_options( argv, "reducer_outfile", "" ) 
if(reducer_outfile == ""): error_exit_with_message("reducer_outfile == \"\"!")
if(dirname(reducer_outfile) != ""):  error_exit_with_message("dirname(reducer_outfile) != \"\"")

delete_files=True

reducer_outfile=data_foldername + '/' + reducer_outfile

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )

deleting_files_signal_file="%s/DAG_rebuild_bulge_reducer_deleting_files_signal_file.txt" %(data_foldername)

if(exists(deleting_files_signal_file)):	error_exit_with_message("%s exist!" %(deleting_files_signal_file))

##################################################################################################################

#######REWROTE THIS ON FEB 08, 2012##########

if( exists(reducer_outfile) ): 
	print "Warning reducer_outfile (%s) already exist....removing!" %(reducer_outfile)
	submit_subprocess("rm %s" %(reducer_outfile))

########################################################

if(exists(condor_submit_file)==False): error_exit_with_message("condor_submit_file (%s) doesn't exist!" %(condor_submit_file) )

lines = safe_readlines(condor_submit_file) 

num_queue_line=0
num_mapper_jobs=0

for line in lines:

	cols=line.split()

	if(len(cols)<2): error_exit_with_message("len(cols)<2 for line (%s)" %(line))

	if(cols[0] != 'Queue'): continue

	num_queue_line+=1

	if( len(cols)!=2 ): error_exit_with_message("len(cols)!=2 for Queue line=(%s)" %(line))

	num_mapper_jobs=int(cols[1])

if(num_queue_line!=1): error_exit_with_message("num_queue_line=(%s)!=1" %(num_queue_line))


########################################################
non_empty_silent_file_list=[]

for queue_ID in range(num_mapper_jobs):

	#FINAL_REBUILD_BULGE/STANDARD/S_8847/mapper_rebuild_bulge.out

	silent_file= "%s/S_%d/mapper_rebuild_bulge.out" %(data_foldername, queue_ID)

	if(exists(silent_file)==False): error_exit_with_message("silent_file (%s) doesn't exist!" %(silent_file))

	silent_data=safe_open(silent_file, mode='r', Is_master=False)

	first_silent_line=silent_data.readline();

	silent_data.close()

	Is_empty=False

	if(first_silent_line=="DAG_builge_bulge.py: Unsucessfully rebuild bulge(s)!!\n"): #No virtaul_ribose

		Is_empty=True

	else:

		assert_is_valid_non_empty_silent_file(silent_file)		

		non_empty_silent_file_list.append(silent_file)

	print "queue_ID=%d | silent_file=%s | Is_empty=%s" %(queue_ID , silent_file, Is_empty)

print "num_mapper_jobs=%d | len(non_empty_silent_file_list)=%d" %(num_mapper_jobs, len(non_empty_silent_file_list))

sys.stdout.flush()
sys.stderr.flush()

if(len(non_empty_silent_file_list)==0): error_exit_with_message("len(non_empty_silent_file_list)==0!")
	
concatenate_outfiles(infile_list=non_empty_silent_file_list, outfile=reducer_outfile, add_file_num_to_tag=False)

if(exists(reducer_outfile)==False): error_exit_with_message("reducer_outfile (%s) doesn't exist!" %(reducer_outfile) )


'''//PRE Feb 08, 2012:
globstring = data_foldername+'/S_*/' + 'mapper_rebuild_bulge.out' 
print "globstring= ", globstring
sys.stdout.flush()

globfiles = glob( globstring )
if(len( globfiles ) == 0):
	sleep( 10 ) 
	globfiles = glob( globstring )
globfiles.sort()

if(len( globfiles ) == 0): error_exit_with_message("ERROR! len(globfiles)==0")


for n in range(len(globfiles)):

	print "globfiles[%d]=%s" %(n, globfiles[n])
'''


##################################################################################################################

if(delete_files):
	create_generic_done_signal_file(deleting_files_signal_file)		#The mark the point where the scipt is now not rerunnable!

	for queue_ID in range(num_mapper_jobs):

		foldername="%s/S_%d/" %(data_foldername, queue_ID)
		if(foldername==""): error_exit_with_message("foldername==\"\"")
		delete_command="rm -r %s " %(foldername)
		print delete_command
		sys.stdout.flush()
		sys.stderr.flush()
		submit_subprocess( delete_command ) 

print_title_text("Exit " + python_command)

