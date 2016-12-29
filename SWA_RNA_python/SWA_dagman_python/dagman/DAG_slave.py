#!/usr/bin/env python

from os import system
from sys import argv
from time import sleep
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options

######################################################################


#######################################################################

START_argv=copy.deepcopy(argv)

print 'Enter: ', list_to_string(START_argv)	#Important to input at least one line into outfile to let master_node know that this slave_job is not broken! Aug 10th, 2011.

sys.stdout.flush()
sys.stderr.flush()

jobdir = parse_options( argv, "slave_dir", "" )

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )

if(jobdir==""): error_exit_with_message("jobdir==\"\"" )

if( exists( jobdir ) == False): error_exit_with_message("jobdir (%s) doesn't exist!" %(done_signal_filename) )

finished_file = jobdir + '/finished.txt'


while(True):

	if(exists(finished_file)):
		print "received finished_file signal (%s). EXITING..."
		break

	command_file_name = jobdir + '/run_this_script.txt'


	if exists( command_file_name ):
		run_first_job=True

		sleep( 5 ) # stupid file locking.....POTENTIAL BUG HERE..assume that after 5 seconds, the file will be written....

		lines = safe_open( command_file_name , mode='r').readlines() 

		if(len(lines)!=2): 
			print "lines= " , lines
			error_exit_with_message("len(lines)!=2")
		
		job_script_filename=lines[0][:-1]
		done_signal_filename=lines[1][:-1]

		if( exists(done_signal_filename)==True ):
			error_exit_with_message("done_signal_filename (%s) already exist!" %(done_signal_filename) )

		
		job_script_lines = safe_open( job_script_filename, mode='r' ).readlines() 

		JS_line= job_script_lines[0][:-1]
		cols = string.split( JS_line )

		print "potential wrapper line_cols= " , cols  

		if (cols[0] ==  "wrapper_script"):
			DAG_reducer_EXE=get_PYEXE("dagman/DAG_reducer.py")

			if( cols[1] != "=" ): error_exit_with_message('cols[0] ==  "wrapper_script" but cols[1] != "="')
			if( cols[2] != DAG_reducer_EXE): error_exit_with_message('cols[0] ==  "wrapper_script" but cols[2] != "%s"' %(DAG_reducer_EXE))
			wrapper_script = cols[2]

			command= wrapper_script + ' ' + job_script_filename
			
		else: #mapper job

			command = 'source '+ job_script_filename

		print command
		submit_subprocess( command )

		print "-----" 
		
		#OK if reach this point means that job was successfully ran. Create the done_signal file;
		fid = safe_open( done_signal_filename, 'w' )
		fid.write( 'MESSAGE FROM SLAVE_NODE at jobdir %s: JOB %s SUCCESSFULLY RAN!\n ' %(jobdir, job_script_filename) )
		fid.close()

	  # After done, remove the script file...
	  # A more robust signal might be to create a "done" file.
		submit_subprocess( 'rm -rf ' + command_file_name )


	sleep( 2 )

print "SUCCESSFULLY RAN: %s" %(list_to_string(START_argv))

