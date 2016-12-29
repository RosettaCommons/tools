#!/usr/bin/env python

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
from os import system
from os.path import basename, dirname, exists, expanduser, abspath
from time import sleep
#import popen2


####COMMONLY CALLED FUNCTIONS################################


def error_exit_with_message(message):
	print >> sys.stdout, "A error had occur, see errfile.txt"
	print >> sys.stderr, "Error!: "  + message

#	traceback.print_exc()
	sys.stdout.flush()
	sys.stderr.flush()

	assert(False)

	sys.exit(1)

def error_check(retcode, command, Is_master=False):

	if retcode is not 0:
		sys.stderr.write("Child was terminated by signal %d \n" %retcode)
		sys.stderr.write("Error subprocess: %s \n" %command)
		sys.stdout.flush()
		sys.stderr.flush()
		if(Is_master):
			kill_all_slave_jobs_and_exit()
		else:
			error_exit_with_message("subprocess submission error")

def submit_subprocess(command, Is_master=False):

	sys.stdout.flush()
	sys.stderr.flush()

#	retcode = subprocess.call(command , shell=True) #This is not avialable in python 2.3.4 (current Biox's python version) might have to make shell false on BIOX??
	retcode = system(command)
#	print "retcode= ", retcode
	error_check(retcode, command, Is_master)
	sys.stdout.flush()
	sys.stderr.flush()

#PIPE on POSIX are implemented using pipe system call so it will have all
#the limitations of the system buffering that is used for implementing
#pipes. The buffer size is system dependent. For example, on my Linux, I
#hit the hang condition at 64K.What basically happens is that the child
#fills up the pipe buffer and is blocked until space is available.

#To prevent memory limitation, I decided to temporary write the data to disk, open() and readlines() and then remove the temp data.
#Compare to strand popen(), this is slower but should be less prone to errors.

#Lastly since the readlines() function is invoked...this function should not be called to read in large files

def	popen_and_readlines(command, Is_master=False):

	sys.stdout.flush()

	if(exists( 'temp_data.txt')):
		submit_subprocess('rm temp_data.txt', Is_master)

	submit_subprocess( command + ' > temp_data.txt', Is_master)

	try:
		lines = open( 'temp_data.txt' ).readlines()
	except:
		if(Is_master):
			kill_all_slave_jobs_and_exit(error_message)
		else:
			error_exit_with_message(error_message)

	submit_subprocess('rm temp_data.txt', Is_master)

	sys.stdout.flush()

	return lines


###Alternative method, faster but prone to memory issues
#	pipe=popen2.Popen3(command)
#	retcode=pipe.wait()
#	error_check(retcode, Is_master)
#	return pipe.fromchild.readlines()

##############################################################################################

def safe_open(filename, mode='r' ,Is_master=True): ###For now assume only master job will call this function

	if(exists( filename ) == False and mode=='r'): ###Sometime Biox file system is slow to respond..so give it sometime
		sleep(20)

	try:
		data = open( filename, mode)
	except:
		error_message="cannot open %s " %filename
		if(Is_master):
			kill_all_slave_jobs_and_exit(error_message)
		else:
			error_exit_with_message(error_message)
	return data

def	line_counts(filename, Is_master=True): ###For now assume only master job will call this function

	data = safe_open(filename, 'r', Is_master)

	count=0
	for line in data:
		count+=1

	data.close()

	return count

##################################################################################################

def print_title_text(title):

	title_length=len(title)
	char_per_line=80;
	dash_length=char_per_line-title_length;

	print

	for i in range( dash_length/2 ):
		sys.stdout.write("-")

	sys.stdout.write(title)

	for i in range( dash_length/2 ):
		sys.stdout.write("-")

	print

	sys.stdout.flush()
	sys.stderr.flush()

##################################################################################################

def get_PYDIR():

	HOMEDIR = expanduser('~')
	PYDIR = HOMEDIR+'/SWA_dagman_python2'
	assert( exists( PYDIR ) )
	return PYDIR
##################################################################################################
def Is_column_name_line(line):

	if (line.find( "description" ) != -1): return True


#def Is_silent_file_data_line(line): #input a string

#	if (line.find( "ANNOTATED_SEQUENCE") != -1): return True
#	if (line.find( "SEQUENCE" ) != -1): return False
#	if (line.find( "description" ) != -1): return False
#	if (line.find( "REMARK BINARY_SILENTFILE" ) != -1): return False

#	return True
##################################################################################################
def get_python_command(argv):
	command=''

	for i in range(len(argv)):
		command+= argv[i] + ' '

	return command

##################################################################################################
def get_silent_file_col_index(col_name, silent_file):

	firstlines = popen_and_readlines('head -n 3 '+ silent_file)

	col_name_line=firstlines[1]
	col_name_list=string.split(col_name_line)

	try:
		col_index=col_name_list.index(col_name)
	except:
		error_exit_with_message("Cannot find col_index with corresponding colname= %s " %(col_name) )

	return col_index
##################################################################################################
def count_struct(silent_file):

	print_title_text("count_struct for silent_file: %s " %(silent_file) )
	col_index=get_silent_file_col_index("score", silent_file)

	SCORE_silentfile='temp.txt'

	command = 'grep SCORE %s > %s' %(silent_file, SCORE_silentfile)

	submit_subprocess(command)

	try:
		infile = open(SCORE_silentfile, 'r' )
	except:
		error_exit_with_message("cannot open %s " % SCORE_silentfile)

	num_struct=0
	first_line=True
	for line in infile:

		if(Is_column_name_line(line)==True): continue ##Take this off it is significantly slow down code

		cols = string.split( line )

		#Check that can read score from score_line
		try:
			score = float( cols[ col_index ] )
		except:
			print cols
			error_exit_with_message("problem with reading score line: %s")
		num_struct+=1

	infile.close()
	submit_subprocess("rm %s" %SCORE_silentfile) #Mod this: Don't remove SCORE_silentfile. Keep it for consistency, doesn't take up too much space anyways.

	return num_struct

def get_min_max(silent_file):


	col_index=get_silent_file_col_index("score", silent_file)
	print "col_index=", col_index
	min_score=0
	max_score=0

	SCORE_silentfile='temp.txt'

	command = 'grep SCORE %s > %s' %(silent_file, SCORE_silentfile)

	submit_subprocess(command)


	try:
		infile = open(SCORE_silentfile, 'r' )
	except:
		error_exit_with_message("cannot open %s " % SCORE_silentfile)

	first_line=True
	for line in infile:

		if(Is_column_name_line(line)==True): continue ##Take this off it is significantly slow down code

		cols = string.split( line )

		try:
			score = float( cols[ col_index ] )
		except:
			print cols
			error_exit_with_message("problem with reading score line: %s")

		if(first_line==True):
			min_score=score
			max_score=score
			first_line=False

		if(min_score > score): min_score=score
		if(max_score < score): max_score=score

	infile.close()
	submit_subprocess("rm %s" %SCORE_silentfile) #Mod this: Don't remove SCORE_silentfile. Keep it for consistency, doesn't take up too much space anyways.

	return (min_score, max_score)

