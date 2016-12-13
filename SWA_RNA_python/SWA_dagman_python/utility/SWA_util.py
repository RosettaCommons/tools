#!/usr/bin/env python


###############################Aug 11th, 2011##################
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
###############################################################

####### SWA reviaval Oct. 2, 2014 #############################
import subprocess

###############################################################

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
from os import system
from os.path import basename, dirname, expanduser, abspath
from os.path import exists as os_path_exists
from time import sleep
import popen2
import copy
from sets import Set
#######################################################


from list_operations import *
from error_util import *
from silent_file_util import *
from PATHS import get_PYDIR, get_PYEXE, get_PYPATH, get_rosetta_EXE, get_rosetta_database_folder
from PATHS import get_rosetta_EXE_specified_no_graphic_string, PATH_exists, is_release_mode, use_new_src_code
#######################################################


####################################################################################

def exists( filename ):
	### This overload is important for preventing false negatives
	attempt = 0
	while attempt < 2:
		status = os_path_exists( filename )
		if status:  break
		attempt += 1
	return status

####################################################################################

def submit_subprocess_allow_retry(command, Is_master=False):

	num_try_max=10

	for n in range(1, num_try_max+1 ):

		sys.stdout.flush()
		sys.stderr.flush()

		if(n==num_try_max):
			submit_subprocess(command, Is_master)
		else:
			#retcode = subprocess.call(command + ' 2> submit_subprocess_retry_err.txt', shell=True)
			retcode = system(command + ' 2> submit_subprocess_retry_err.txt')
			retcode = retcode % 256 # exit values greater than 255 return an exit code modulo 256 (added Feb. 23, 2015)

			if(retcode==0):
				break
			else:
				print "retcode!=0 (%d) for command %s !!" %(retcode, command)
				fid = open( 'PROBLEMATIC_RETCODE_SUBPROCESS.LOG', 'a' )
				fid.write( "retcode!=0 (%d) for command %s !! \n" %(retcode, command) )
				fid.close()

		sleep(2)
		sys.stdout.flush()
		sys.stderr.flush()

	sys.stdout.flush()
	sys.stderr.flush()


#######################################################

def submit_subprocess(command, Is_master=False):

	sys.stdout.flush()
	sys.stderr.flush()

	#retcode = subprocess.call(command , shell=True) #This is not avialable in python 2.3.4 (current Biox's python version) might have to make shell false on BIOX??
	retcode = system(command)
	retcode = retcode % 256 # exit values greater than 255 return an exit code modulo 256 (added Feb. 23, 2015)


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

#######################################################

def safe_readlines(filename): #June 13, 2011.

	if(exists(filename)==False): error_exit_with_message("filename (%s) doesn't exist!" %(filename) )

	try:
		line_list = open( filename  ).readlines()
	except:
		error_message="cannot open and readlines filename (%s) " %(filename)
		error_exit_with_message(error_message)

	return line_list

#######################################################


def	popen_and_readlines(command, Is_master=False, tag=""):

	sys.stdout.flush()

	if(Is_master): #problem is now both master and post_process slave call the popen_and_readlines function!
		if(tag!=""): error_exit_with_message('Is_master=True but tag!=""!' )
		temp_data_filename='master_temp_data.txt'
	else: #called by slave
		if(tag==""): error_exit_with_message('Is_master=False but tag==""!' )
		if(tag.count('/')!=0): error_exit_with_message('Invalid tag %s since contain "/"!' %(tag) )
		temp_data_filename='%s_temp_data.txt' %(tag)
		#print "temp_data_filename= %s " %(temp_data_filename)


	if(exists( temp_data_filename)):
		sleep(5)
		if(exists( temp_data_filename)): error_exit_with_message("temp_data_filename %s already exists!" %(temp_data_filename) )

	submit_subprocess_allow_retry( command + ' > %s' %(temp_data_filename), Is_master)

	lines=safe_open(temp_data_filename, 'r' ,Is_master).readlines()   #Sept 30, 2010

	if(exists( temp_data_filename)==False): error_exit_with_message("temp_data_filename %s does not exist!" %(temp_data_filename) )
	submit_subprocess_allow_retry('rm %s' %(temp_data_filename), Is_master)

	sys.stdout.flush()

	return lines


###Alternative method, faster but prone to memory issues
#	pipe=popen2.Popen3(command)
#	retcode=pipe.wait()
#	error_check(retcode, Is_master)
#	return pipe.fromchild.readlines()

##############################################################################################
'''
def safe_open(filename, mode='r' ,Is_master=True):

	directory_name=dirname( filename )
	if( (directory_name !="") and (exists( directory_name ) == False) ):
		print "creating local directory:  %s" %(directory_name)
		submit_subprocess_allow_retry('mkdir -p %s ' %(directory_name), Is_master)

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
'''


def safe_open(filename, mode='r' ,Is_master=True):

	if(mode!='r'):
		directory_name=dirname( filename )
		if( (directory_name !="") and (exists( directory_name ) == False) ):
			print "creating local directory:  %s" %(directory_name)
			submit_subprocess_allow_retry('mkdir -p %s ' %(directory_name), Is_master)

	if(exists( filename ) == False and mode=='r'): ###Sometime Biox file system is slow to respond..so give it sometime
		sleep(20)

	for n in range(10):
		try:
			data = open( filename, mode)
			return data
		except:
			print "problem opening file %s , mode=%s .....RETRYING...." %(filename, mode)
			sleep(2) #March 16, 2011

	#if reach this point, means that code have failed..
	error_message="cannot open %s mode=%s" %(filename, mode)
	if(Is_master):
		kill_all_slave_jobs_and_exit(error_message)
	else:
		error_exit_with_message(error_message)



##############################################################################################

def local_safe_mkdir(foldername, Is_master=True):

	if( foldername!="" and exists(foldername)==False): submit_subprocess_allow_retry('mkdir -p ' + foldername, Is_master)

##############################################################################################

def	line_counts(filename, Is_master=True):

	data = safe_open(filename, 'r', Is_master)

	count=0
	for line in data:
		count+=1

	data.close()

	return count

##################################################################################################

def print_title_text(title):

	title_length=len(title)
	char_per_line=207;
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
def Is_column_name_line(line):

	if (line.find( "description" ) != -1): return True

	return False # add this line on April 10th, 2011

#def Is_silent_file_data_line(line): #input a string

#	if (line.find( "ANNOTATED_SEQUENCE") != -1): return True
#	if (line.find( "SEQUENCE" ) != -1): return False
#	if (line.find( "description" ) != -1): return False
#	if (line.find( "REMARK BINARY_SILENTFILE" ) != -1): return False

#	return True
##################################################################################################
def get_python_command(argv_local):
	command=''

	for i in range(len(argv_local)):
		command+= argv_local[i] + ' '

	return command


##################################################################################################
def get_silent_file_col_index(col_name, silent_file):

	firstlines = popen_and_readlines('head -n 3 '+ silent_file, Is_master=False, tag="local")

	try:
		col_name_list=firstlines[1].split() # standard silent_file
		col_index=col_name_list.index(col_name)
	except:
		try:
			col_name_list=firstlines[0].split() #grep SCORE: file
			col_index=col_name_list.index(col_name)
		except:
			error_exit_with_message("Cannot find col_index with corresponding colname= %s " %(col_name) )

	return col_index


################################################################
def get_min_max(silent_file):


	col_index=get_silent_file_col_index("score", silent_file)
#	print "col_index=", col_index
	min_score=0
	max_score=0

	SCORE_silentfile='temp.txt'

	command = 'grep "SCORE: " %s > %s' %(silent_file, SCORE_silentfile)

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


#####Used by DAG_LFS_continuous.py and DAG_reducer.py
def	get_job_info(log_foldername, job_name, CreateDirs = True ):

	job_info = {}

	prefix=""
	if(job_name=='reducer' or job_name=='0'):
		prefix='KEEP_LOG_FILE/'

	job_info['out_log_filename'] = prefix + log_foldername + '/outfile/' + job_name + '.out'
	job_info['err_log_filename'] = prefix + log_foldername + '/errfile/' + job_name + '.err'
	job_info['job_script_filename'] = prefix + log_foldername + '/job_script/' + job_name + '.qsub'
	job_info['done_signal_filename'] = 'DONE/' + log_foldername + '/' + job_name + '.txt'

	if CreateDirs:
		for filename in  [ 'out_log_filename', 'err_log_filename', 'job_script_filename', 'done_signal_filename']:
			outdir = dirname( job_info[ filename] )
			if not exists(  outdir ): submit_subprocess_allow_retry('mkdir -p %s ' %(outdir), True )


	return job_info

#####Used by DAG_LFS_continuous.py and DAG_reducer.py


def check_if_job_in_done( job_info_list ):

	if(len(job_info_list)==0):
		print "Jobs is done: len(job_info_list)==0)"
		return True

	num_still_running_job=0

	for job_info in job_info_list:

		done_signal_filename=job_info['done_signal_filename']
#		print "done_singal_filename= ", done_signal_filename
#		if( (exists(done_signal_filename)== False) and (exists(done_signal_filename+'\n')== False) ) : num_still_running_job +=1
		if( (exists(done_signal_filename)== False) ) : num_still_running_job +=1

	print_job_string=""
	for n in range(len(job_info_list)-1, -1, -1): #goes from len(job_info_list)-1 to 0
		print_job_string=job_info_list[n]['done_signal_filename']
		if(print_job_string.count('reducer')==0): break


	if(num_still_running_job==0):
		print "Jobs is done: " , print_job_string
		return True
	else:

		print "Jobs still running: ", print_job_string, ' ', num_still_running_job
		return False



####################################################################
def check_valid_VDW_rep_screen_info_list(VDW_rep_screen_info_list):

	if( len(VDW_rep_screen_info_list)== 0): error_exit_with_message("len(VDW_rep_screen_info)== 0" )
	if( len(VDW_rep_screen_info_list[0])== 0): return

    ################################################################
	### VDW_rep_screen_info_list no longer requires additional arguments,
	### user only needs to provide a pdb or a list of pdbs
	###
	if( len( VDW_rep_screen_info_list ) == 1 ): return True
	if( len( VDW_rep_screen_info_list )  > 1 ):
		### Check to see if VDW_rep_screen_info_list is a list of pdbs.
		if( sum([ ('.pdb' in info) for info in VDW_rep_screen_info_list ]) == len( VDW_rep_screen_info_list ) ): return True
    ###
 	### -- caleb, 11.19.2014
    ################################################################

	if( (len(VDW_rep_screen_info_list) % 3) != 0):
		print "VDW_rep_screen_info_list: ", VDW_rep_screen_info_list
		error_exit_with_message("len(VDW_rep_screen_info_list) % 3 != 0")

	for n in range(len(VDW_rep_screen_info_list)):
		if(n % 3 == 0):
			if(exists( VDW_rep_screen_info_list[n] )==False): error_exit_with_message("'exists( VDW_rep_screen_info_list[%d] (%s) )==False" %(n, VDW_rep_screen_info_list[n]) )


####################################################################
def Is_valid_empty_silent_file(silent_file, verbose=True):

	prefix_reason_string="silent_file (%s) is not a valid_empty_silent_file." %(silent_file)

	if ( not exists(silent_file) ):
		if (verbose): 
			print "%s REASON: silent_file doesn't exist!" % (prefix_reason_string)
		return False

	data = safe_open(silent_file, 'r', Is_master=False)
	data_lines = data.readlines()
	data.close()

	if ( len(data_lines) < 1 ):
		if (verbose): 
			print "%s REASON: empty silent_file num_lines < 1!" % (prefix_reason_string)
		return False

	first_line = data_lines[0]
	if (verbose):
			print "first_line= ", first_line

	### EXCEPTIONS HERE
	valid_exceptions = [
		"empty cluster silent_file since all input_silent_file are empty.",
		"empty filtered silent_file since no non-empty sampler silent_file.",
		"Empty filterer_outfile. No struct_pair passed screen.",
		"empty cluster silent_file since at least one of the two input_silent_file is empty.",
		"empty cluster silent_file since at least of the two input_silent_file is empty."
	]

	return ( any ( exception in first_line for exception in valid_exceptions ) )
		

####################################################################
def Is_valid_non_empty_silent_file(silent_file, verbose=True):

	prefix_reason_string="silent_file (%s) is not a valid_non_empty_silent_file." %(silent_file)

	if(exists(silent_file)==False):
		if(verbose): print "%s REASON: silent_file doesn't exist!" %(prefix_reason_string)
		return False

	####OK FIRST CHECK THAT THERE ARE AT LEAST THREE LINES IN SILENT_FILE#####
	data = safe_open(silent_file, 'r', Is_master=False)

	line_num=0

	SEQUENCE_LINE=""
	COLUMN_NAME_LINE=""

	for line in data:

		line_num+=1

		if(line_num==1): SEQUENCE_LINE=line # The SEQUENCE: gggcgcagccu line

		if(line_num==2): COLUMN_NAME_LINE=line # The column name line

		if(line_num==3): break

	data.close()

	if(line_num!=3):
		if(verbose): print "%s REASON: num_lines=(%s)<3!" %(prefix_reason_string, line_num)
		return False

	COL_NAME_LIST=COLUMN_NAME_LINE.split()

	assert_no_duplicate_in_string_list(COL_NAME_LIST)

	if(SEQUENCE_LINE[0:9]!='SEQUENCE:'):
		if(verbose): print "%s REASON: SEQUENCE_LINE[0:9]!='SEQUENCE:' for SEQUENCE_LINE (%s)" %(prefix_reason_string, SEQUENCE_LINE)
		return False

	if(COLUMN_NAME_LINE[0:6]!='SCORE:'):
		if(verbose): print "%s REASON: COLUMN_NAME_LINE[0:6]!='SCORE:'" %(prefix_reason_string)
		return False

	if(COLUMN_NAME_LINE.count('description') != 1 ):
		if(verbose): print "%s REASON: COLUMN_NAME_LINE.count('description') != 1" %(prefix_reason_string)
		return False

	if(COL_NAME_LIST[0]!='SCORE:'):
		if(verbose): print "%s REASON: COL_NAME_LIST[0]!='SCORE:'" %(prefix_reason_string)
		return False


	if(COL_NAME_LIST[-1]!='description'):
		if(verbose): print "%s REASON: COL_NAME_LIST[-1]!='description'" %(prefix_reason_string)
		return False

	return True

####################################################################

def assert_is_valid_non_empty_silent_file(silent_file, verbose=True):

	if( Is_valid_non_empty_silent_file(silent_file, verbose)==False):
		error_exit_with_message("silent_file (%s) is not a valid_non_empty_silent_file!" %(silent_file))

####################################################################

##Feb 08, 2012:THIS FUNCTION IS DEPRECATE, consider using Is_valid_non_empty_silent_file or assert_is_valid_non_empty_silent_file() instead!
def is_non_empty_silent_file(silent_file, verbose=False):

	if(exists(silent_file)==False): return False

	data = safe_open(silent_file, 'r', Is_master=False)

	count=0
	for line in data:
		count+=1
		if(count==3): break

	data.close()

	if(count==3):
		return True
	else:
		if(verbose): print "silent_file (%s) exists but is empty!" %(silent_file)
		return False


##################################################################################################
def count_struct(silent_file):

	#print_title_text("count_struct for silent_file: %s " %(silent_file) )

	assert_is_valid_non_empty_silent_file(silent_file)

	col_index=get_silent_file_col_index("score", silent_file)

	SCORE_silentfile='temp.txt'

	command = 'grep "SCORE: " %s > %s' %(silent_file, SCORE_silentfile)

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
	submit_subprocess("rm %s" %(SCORE_silentfile))

	return num_struct
####################################################################


def append_to_basename(suffix, filename):

	move_out_one_folder=False

	if( (len(suffix)>3) and (suffix[0:3]=="../") ):
		move_out_one_folder=True
		suffix=suffix[3:]

	NEW_filename=suffix + basename(filename)

	if(dirname(filename)!=""):
		directory_name=dirname(filename)

		if(move_out_one_folder):
			right_slash_index=directory_name.rfind("/")
			if(right_slash_index==-1): error_exit_with_message("right_slash_index==-1")
			directory_name=directory_name[:right_slash_index]

		NEW_filename=directory_name + "/" + NEW_filename
	else:
		if(move_out_one_folder): error_exit_with_message("move_out_one_folder==True but dirname(filename)==\"\"!")

	return NEW_filename


####################################################################

def Is_valid_value_dinucleotide_at_single_element_cc(dinucleotide_at_single_element_cc):

	if(dinucleotide_at_single_element_cc!="all" and dinucleotide_at_single_element_cc!="pure_append_prepend" and dinucleotide_at_single_element_cc!="none"):
		error_exit_with_message("Invalid dinucleotide_at_single_element_cc value=%s" %(dinucleotide_at_single_element_cc) )


####################################################################

def sort_dictionary_list(dict_list, key): #dict_list is actually a reference to the actual dict_list object!

	######Python 2.4 and after way!####
	#dict_list=sorted(dict_list, key=lambda dict_object: dict_object[key] )

	######Prepython 2.4 way!####
	decorated = [(dict_object[key], i, dict_object) for (i, dict_object) in enumerate(dict_list)]

	decorated.sort()

	dict_list=[tuple_object[2] for tuple_object in decorated] #OH, here the dict_list reference now point to a new object!

	decorated=[] #release memory??? #But it just pointing to a new object??

	return dict_list #need to return new dict_list reference since it now point to a different object!

####################################################################
def ensure_no_duplicate_tags(silent_file):

	print "Ensuring no_duplicate_tags in %s...." %(silent_file)

	data=safe_open(silent_file, mode='r', Is_master=False)

	SEQUENCE_LINE=data.readline()
	COLUMN_NAME_LINE=data.readline()

	COL_NAME_LIST=COLUMN_NAME_LINE.split()

	try:
		tag_col_index=COL_NAME_LIST.index('description')
	except:
		print "COL_NAME_LIST=", COL_NAME_LIST
		error_exit_with_message("Cannot find description column index!")

	if(tag_col_index!=(len(COL_NAME_LIST)-1)): error_exit_with_message("tag_col_index!=(len(COL_NAME_LIST)-1)")

	tag_map={}

	while(True):

		line=data.readline()

		if(line==''): break #End of file!

		if(len(line) <= 1): error_exit_with_message("len(line) <= 1") #check for line with only '\n'

		if(line.count('SCORE:') == 0): continue

		if(line[0:6]!="SCORE:"): error_exit_with_message("line.count('SCORE:') != 0 but line[0:6]!=\"SCORE:\" for line=%s" %(line))

		if(line.find('description') != -1): error_exit_with_message("extra column_name line (%s)" %(line) )

		line_split=line.split()

		if(tag_col_index!=(len(line_split)-1)): error_exit_with_message("tag_col_index!=(len(line_split)-1)) for line=%s" %(line))

		tag=line_split[tag_col_index]

		if(tag_map.has_key(tag)): error_exit_with_message("A prior score_line already used the tag (%s)!" %(tag) )

		tag_map[tag]=True #Added this on Dec 20, 2011! Previously, missing this crucial line

	data.close()


#####################################################
def convert_string_to_float(string_num):

	Is_negative=False

	if(string_num[0]=="N"):
		#print "string_num= ", string_num
		string_num=string_num[1:]
		Is_negative=True

	float_num=999.99

	try:
		float_num=float(string_num)
	except:
		error_exit_with_message("Unable to parse the string_num (%s)" %(string_num))

	if(Is_negative):
		return (-1.0*float_num)
	else:
		return ( 1.0*float_num)

#####################################################
def create_generic_done_signal_file(file_name):

	if(exists(file_name)==True): error_exit_with_message("done_signal_file_name (%s) already exist!" %(file_name))

	outfile = safe_open(file_name, mode='w' ,Is_master=False)

	outfile.write("GENERIC_DONE_SIGNAL\n" )
	outfile.close()

##############################################
def check_tag( check_string, tag ): # allow hash tags
	return ( check_string == tag or check_string == "#"+tag )

