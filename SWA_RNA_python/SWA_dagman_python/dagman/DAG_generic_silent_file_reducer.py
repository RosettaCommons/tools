#!/usr/bin/env python


from os import system,popen
from os.path import dirname,basename,expanduser,abspath
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


#~/SWA_RNA_python/SWA_dagman_python/dagman/DAG_generic_silent_file_reducer.py -reducer_outfile minimize_WITH_SHIFT_STATS_Jan_23_SWA.out -outfolder FULL_LENGTH_MINIMIZER -condor_submit_file CONDOR//TEST_SWA_MINIMIZE.condor   -delete_mapper_files False

#~/SWA_RNA_python/SWA_dagman_python/dagman/DAG_generic_silent_file_reducer.py -reducer_outfile minimize_WITH_SHIFT_STATS_Jan_23_SWA.out -outfolder FULL_LENGTH_MINIMIZER -condor_submit_file CONDOR//SWA_MINIMIZE.condor  

python_command = get_python_command(argv)
print_title_text("Enter: " + python_command)

data_foldername= parse_options( argv, "outfolder", "" ) 

if(data_foldername==""): error_exit_with_message("User need to pass in data_foldername ")

if(exists( data_foldername )==False): error_exit_with_message("data_foldername (%s) doesn't exist!" %(abspath(data_foldername) ))

condor_submit_file=parse_options( argv, "condor_submit_file", "" ) #FEB 08, 2012
if(condor_submit_file == ""): error_exit_with_message("condor_submit_file == \"\" ") 

reducer_outfile= parse_options( argv, "reducer_outfile", "" ) 
if(reducer_outfile == ""): error_exit_with_message("reducer_outfile == \"\"!")
if(dirname(reducer_outfile) != ""):  error_exit_with_message("dirname(reducer_outfile) != \"\"")

delete_mapper_files= parse_options( argv, "delete_mapper_files", "True" ) 

reducer_outfile=data_foldername + '/' + reducer_outfile

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )

deleting_files_signal_file="%s/DAG_rebuild_bulge_reducer_deleting_files_signal_file.txt" %(data_foldername)

if(exists(deleting_files_signal_file)):	error_exit_with_message("%s exist!" %(deleting_files_signal_file))

##################################################################################################################

if( exists(reducer_outfile) ): 
	print "Warning reducer_outfile (%s) already exist....removing!" %(reducer_outfile)
	submit_subprocess("rm %s" %(reducer_outfile))

########################################################

if(exists(condor_submit_file)==False): error_exit_with_message("condor_submit_file (%s) doesn't exist!" %(condor_submit_file) )

lines = safe_readlines(condor_submit_file) 

num_queue_line=0
num_mapper_outfiles_line=0

num_mapper_jobs=0
mapper_outfile=""

for line in lines:

	cols=line.split()

	if(len(cols)<2): error_exit_with_message("len(cols)<2 for line (%s)" %(line))

	if(cols[0] == 'Queue'): 

		num_queue_line+=1

		if( len(cols)!=2 ): error_exit_with_message("len(cols)!=2 for Queue line=(%s)" %(line))

		num_mapper_jobs=int(cols[1])

	if( cols[0] == "mapper_outfiles"  ): 

		num_mapper_outfiles_line=+1

		####SOME POSSIBLE mapper_outfiles examples######
		#mapper_outfiles =  FULL_LENGTH_MINIMIZER/$(Process)//minimize_silent_file.out 

		#mapper_outfiles = FINAL_REBUILD_BULGE/STANDARD//S_$(Process)//mapper_rebuild_bulge.out FOLDER:FINAL_REBUILD_BULGE/STANDARD//S_$(Process)/ 

		#For general case, can have more than one mapper_outfile (e.g. second example above...still have to rewrite code to account for this possibility!)
		if( len(cols)!=3 ): error_exit_with_message("len(cols)!=3 for mapper_outfiles line=(%s)" %(line)) 

		if( cols[1] != "=" ): error_exit_with_message("cols[1] != \"=\" for mapper_outfiles line (%s)" %(line))

		mapper_outfile  = cols[2]


if(num_queue_line!=1): error_exit_with_message("num_queue_line=(%s)!=1" %(num_queue_line))
if(num_mapper_outfiles_line!=1): error_exit_with_message("num_mapper_outfiles_line=(%s)!=1" %(num_mapper_outfiles_line))

if(mapper_outfile==""): error_exit_with_message("mapper_outfile=\"\"")
if(num_queue_line<=0): error_exit_with_message("num_queue_line<=0")


########################################################
non_empty_silent_file_list=[]

possible_empty_silent_file_message=[]

possible_empty_silent_file_message.append("DAG_builge_bulge.py: Unsucessfully rebuild bulge(s)!!\n") #No virtaul_ribose
possible_empty_silent_file_message.append("StepWiseRNA_Minimizer:: num_pose_outputted==0, empty silent_file!\n") 
possible_empty_silent_file_message.append("StepWiseMinimizer:: num_pose_outputted==0, empty silent_file!\n") 
##ADD MORE MESSAGE AS NEEDED##

for queue_ID in range(num_mapper_jobs):

	silent_file=mapper_outfile.replace('$(Process)','%d' %(queue_ID) )

	if(exists(silent_file)==False): error_exit_with_message("silent_file (%s) doesn't exist!" %(silent_file))

	silent_data=safe_open(silent_file, mode='r')

	first_silent_line=silent_data.readline();

	silent_data.close()

	Is_empty=False

	if(first_silent_line in possible_empty_silent_file_message):

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

##################################################################################################################
if(delete_mapper_files):
	create_generic_done_signal_file(deleting_files_signal_file)		#The mark the point where the scipt is now not rerunnable [DUE TO delete_mapper_files]!

##Delete mapper_files since they are no longer needed.
for queue_ID in range(num_mapper_jobs):

	silent_file=mapper_outfile.replace('$(Process)','%d' %(queue_ID) )

	foldername=dirname(silent_file)

	if(foldername==""): error_exit_with_message("foldername==\"\"")

	if(delete_mapper_files):
		delete_command="rm -r %s " %(foldername)
		print delete_command
		sys.stdout.flush()
		sys.stderr.flush()
		submit_subprocess( delete_command ) 

print_title_text("Exit " + python_command)

