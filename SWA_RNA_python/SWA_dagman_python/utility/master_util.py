#!/usr/bin/env python

######################################################################

from SWA_dagman_python.utility.SWA_util import *

from SWA_dagman_python.scheduler.scheduler_util import *
######################################################################
#This python script contains the master version of:
#./error_util.py:def error_check(retcode, command):
#./subprocess_util.py:def submit_subprocess_allow_retry(command): 
#./subprocess_util.py:def submit_subprocess(command):
#./SWA_util.py:def safe_mkdir(foldername):
#./SWA_util.py:def       line_counts(filename):
#./io_util.py:def        popen_and_readlines(command, tag):
#./io_util.py:def safe_open(filename, mode='r'): 


######################################################################
def master_error_check(retcode, command):

	if(retcode != 0):
		sys.stderr.write("Child was terminated by signal %d \n" %retcode)
		sys.stderr.write("Error subprocess: %s \n" %command)
		sys.stdout.flush()
		sys.stderr.flush()

		master_kill_all_slave_jobs_and_exit()		

######################################################################

def master_submit_subprocess(command):
	
	sys.stdout.flush()
	sys.stderr.flush()
		
	#retcode = system(command) #Comment out on June 05, 2012
	retcode = subprocess.call(command, shell=True) #Change to this on June 05, 2012

	master_error_check(retcode, command)	
	sys.stdout.flush()
	sys.stderr.flush()

######################################################################
def master_submit_subprocess_allow_retry(command):

	num_try_max=10

	for n in range(1, num_try_max+1 ):

		sys.stdout.flush()
		sys.stderr.flush()

		if(n==num_try_max):
			master_submit_subprocess(command)
		else:
			#retcode = system(command + ' 2> submit_subprocess_retry_err.txt') #Comment out on June 05, 2012
			retcode = subprocess.call(command + ' 2> submit_subprocess_retry_err.txt', shell=True) #Change to this on June 05, 2012

			if(retcode==0):
				break
			else:
				print "retcode!=0 (%d) for command (%s)!! Wait 20 seconds before resubmitting command..." %(retcode, command)

		sleep(20) #Change from 5 to 20 seconds on May 14, 2012 #Change from 2 seconds to 5 seconds on May 11, 2012
		sys.stdout.flush()
		sys.stderr.flush()

	sys.stdout.flush()
	sys.stderr.flush()


######################################################################

def master_safe_mkdir(foldername):

	if( foldername!="" and exists(foldername)==False): master_submit_subprocess_allow_retry('mkdir -p ' + foldername)

######################################################################

def master_safe_open(filename, mode='r' ): 


	if(mode!='r'):
		directory_name=dirname( filename )
		if( (directory_name !="") and (exists( directory_name ) == False) ): 
			print "creating local directory:  %s" %(directory_name)
			master_submit_subprocess_allow_retry('mkdir -p %s ' %(directory_name))

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

	master_kill_all_slave_jobs_and_exit(error_message)

#######################################################

def	master_popen_and_readlines(command):

	sys.stdout.flush()
	
	temp_data_filename='master_temp_data.txt'  

	if(exists(temp_data_filename)):
		sleep(5)
		if(exists( temp_data_filename)): master_kill_all_slave_jobs_and_exit("temp_data_filename (%s) already exist!" %(temp_data_filename) )
		
	master_submit_subprocess_allow_retry( command + ' > %s' %(temp_data_filename))

	lines=master_safe_open(temp_data_filename, 'r').readlines()   #Sept 30, 2010

	if(exists( temp_data_filename)==False): master_kill_all_slave_jobs_and_exit("temp_data_filename (%s) doesn't not exist!" %(temp_data_filename) )
	master_submit_subprocess_allow_retry('rm %s' %(temp_data_filename))		

	sys.stdout.flush()
		
	return lines

######################################################################

def master_line_counts(filename):

	data = master_safe_open(filename, 'r')

	count=0
	for line in data:
		count+=1	

	data.close()

	return count




