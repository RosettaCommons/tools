#!/usr/bin/env python

from os import system,chdir
from os.path import dirname
from sys import argv,exit
from time import sleep
import os
######################################################################


from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
from SWA_dagman_python.dagman.DAG_general_util import get_job_info, check_if_job_in_done

######################################################################

wrapper_command = list_to_string(argv)
print 'Enter: ', wrapper_command


if(len(argv)<2): error_exit_with_message("len(argv)<2")

job_script_file = argv[ 1 ]

lines = safe_open(job_script_file, mode='r').readlines()

command=""
reducer_job_name=""
log_foldername=""
num_mapper_jobs=-1

for line in lines:
	cols = string.split( line )

	if(len(cols)<2): error_exit_with_message("len(cols)<2 for line (%s)" %(line))

	if( cols[0] == "num_mapper_jobs") :
		if( len( cols )!= 2 ): error_exit_with_message("len( cols )!=2 for num_mapper_jobs line (%s)!" %(line))
		num_mapper_jobs = int( cols[1] )
		continue

	####OK, all the other lines should have (cols[1] == "=")####
	if( cols[1] != "=" ): error_exit_with_message("cols[1] != \"=\" for line (%s)" %(line))

	if( cols[0] ==  "command" ): command = string.join(cols[2:])

	if( cols[0] ==  "job_name" ):
		if( len(cols) != 3 ): error_exit_with_message("len(cols) != 3 for line (%s)" %(line))
		reducer_job_name = cols[2]

	if( cols[0] ==  "log_foldername" ):
		if( len(cols) != 3 ): error_exit_with_message("len(cols) != 3 for line (%s)" %(line))
		log_foldername = cols[2]

print "-----------------------------------------------------------------------------------------"
print "command=%s" %(command)
print "job_name=%s" %(reducer_job_name)
print "log_foldername=%s" %(log_foldername)
print "num_mapper_jobs=%d" %(num_mapper_jobs)
print "-----------------------------------------------------------------------------------------"

if(command=="" or reducer_job_name=="" or log_foldername=="" or num_mapper_jobs<0):
	exit_message= "Problem parsing job_script. command=%s job_name=%s log_foldername=%s num_mapper_jobs=%d " %(command, reducer_job_name, log_foldername, num_mapper_jobs)
	error_exit_with_message(exit_message)


####################################################################
reducer_job_info = get_job_info(log_foldername, reducer_job_name)

out_log_filename     =reducer_job_info['out_log_filename']
err_log_filename     =reducer_job_info['err_log_filename']
job_script_filename  =reducer_job_info['job_script_filename'] #this is the file which contain the log_foldername and the reducer job_name
done_signal_filename =reducer_job_info['done_signal_filename'] #creating the done signal is taken care for by the DAG_slave for now..

outfile= safe_open(out_log_filename, mode="w")
errfile= safe_open(err_log_filename, mode="w")

sys.stdout.flush()

os.dup2(outfile.fileno(), sys.stdout.fileno()) #this point sys.stdout to outfile

sys.stderr.flush()

os.dup2(errfile.fileno(), sys.stderr.fileno()) #this point sys.stderr to errfile

print "-----------------------------------------------------------------------------------------"
print "command=%s" %(command)
print "job_name=%s" %(reducer_job_name)
print "log_foldername=%s" %(log_foldername)
print "num_mapper_jobs=%d" %(num_mapper_jobs)
print "-----------------------------------------------------------------------------------------"

##############################Check that all the mapper jobs are finish#################################

#if (num_mapper_jobs == 0):
#	print "EARLY EXIT: num_mapper_jobs=%d " % (num_mapper_jobs)
#	print_title_text("Exit: " + wrapper_command) #Early exit
#	exit(0)

count=0
any_mapper_still_running=True
while any_mapper_still_running:
	count+=1

	print "waiting for mapper jobs to finish, count= %d" %(count)
	sys.stdout.flush()
	sys.stderr.flush()

	#NEED to prevent the reducer from being stuck in a infinite while loop if a error occur at the master node or another slave node...
  ###OK IN LSF mode, this is taken care of by bkill!

	any_mapper_still_running=False

	mapper_job_info_list=[]

	for q in range(num_mapper_jobs):

		mapper_job_name='%d' % q

		mapper_job_info = get_job_info(log_foldername, mapper_job_name)

		mapper_job_info_list.append(mapper_job_info)

	job_is_done = check_if_job_in_done( mapper_job_info_list )

	any_mapper_still_running= (job_is_done==False)

	sleep( 2 )

###########################################################################################
print "-----------------------------------------------------------------------------------------"
print "-----------------------------------------------------------------------------------------"
print "------------------All mapper jobs are done! Submitting reducer command!------------------"
print "-----------------------------------------------------------------------------------------"
print "-----------------------------------------------------------------------------------------"


#command_act= "%s > %s 2> %s\n\n" % ( command, out_log_filename, err_log_filename)
command_act= "%s \n\n" % ( command )


print command_act
submit_subprocess( command_act ) #post_process command
submit_subprocess("echo 'Exit reducer job, please check errfile to ensure that there are no errors.' >> %s \n" % (out_log_filename) )



