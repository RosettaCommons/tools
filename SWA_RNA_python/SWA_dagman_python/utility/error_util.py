#!/usr/bin/env python

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
from os import system
from os.path import basename, dirname, exists, expanduser, abspath
from time import sleep
import popen2
import copy
from sets import Set

##############################################

def error_exit_with_message(message):
	print >> sys.stdout, "A error had occur, see errfile.txt"
	print >> sys.stderr, "Error!: "  + message
	print >> sys.stderr, "#####################################################################"
	print >> sys.stderr, "#####################################################################"
	print >> sys.stderr, "#####################################################################"
	sys.stdout.flush()
	sys.stderr.flush()
	
	assert(False)


##############################################

def kill_all_slave_jobs_and_exit(exit_message=""):
	
	
	JOBDIR = 'SLAVE_JOBS/'
	job_tag = abspath( JOBDIR ).replace('/','_')

	sys.stdout.flush()
	sys.stderr.flush()

	command = 'bkill -J %s*' %job_tag	
	print( command )
	system( command )

#	sleep(10)

	sys.stdout.flush()
	sys.stderr.flush()
	exit_message2="Killed all slave jobs and about to exit master script, " + exit_message
	error_exit_with_message(exit_message2)
	
	
def error_check(retcode, command, Is_master=False):

	if( retcode < 0 or retcode > 255 ):
		retcode = retcode % 256 # exit values greater than 255 return an exit code modulo 256 (added Feb. 23, 2015)

	if(retcode != 0):
		sys.stderr.write("Child was terminated by signal %d \n" %retcode)
		sys.stderr.write("Error subprocess: %s \n" %command)
		sys.stdout.flush()
		sys.stderr.flush()
		if(Is_master):
			kill_all_slave_jobs_and_exit()		
		else:
			error_exit_with_message("subprocess submission error")
