#!/usr/bin/env python

from sys import argv,exit,stdout
import string
from glob import glob
from os import system,popen,getcwd
from os.path import basename,exists,expanduser,dirname,abspath
from time import sleep
from SWA_util import *

#######################################
# This could actually be a class!
#######################################

###########When is this called#######when SWA_dagman_LSF_continuous is imported?###
PYDIR = get_PYDIR()
assert( exists( PYDIR ) )
SLAVE_EXE = PYDIR + '/SWA_dagman_slave.py'
####################################################################################

def get_queued_jobs_status():
	'''qstat -f:
	Job Id: 11982.frontend1.local
		  Job_Name = LONG_LONG_LONG_LONG_AAAAAAAAAAAA_LONG_LONG_LONG_LONG_BBBBBBBBBB
		      B_NAME
		  Job_Owner = sripakpa@frontend1.local
		  job_state = R
		  queue = SP
		  server = frontend1.local
		  Checkpoint = u
		  ctime = Fri Mar 30 22:35:26 2012
		  Error_Path = frontend1.local:/home/sripakpa/QSUB_TEST/parin_test_job.err
		  exec_host = node-2-1.local/7
		  Hold_Types = n
		  Join_Path = n
		  Keep_Files = n
		  Mail_Points = be
		  mtime = Fri Mar 30 22:35:31 2012
		  Output_Path = frontend1.local:/home/sripakpa/QSUB_TEST/parin_test_job.out
		  Priority = 0
		  qtime = Fri Mar 30 22:35:26 2012
	'''

	'''qstat -f:
	#Umm there seem to be some variability in the order of the entries.
	Job Id: 11979.frontend1.local
		  Job_Name = qsub_script_create_file_test.txt
		  Job_Owner = sripakpa@frontend1.local
		  resources_used.cput = 00:00:00
		  resources_used.mem = 0kb
		  resources_used.vmem = 0kb
		  resources_used.walltime = 00:00:00
		  job_state = E
		  queue = SP
		  server = frontend1.local
		  Checkpoint = u
		  ctime = Fri Mar 30 21:24:36 2012
		  Error_Path = frontend1.local:/home/sripakpa/QSUB_TEST/qsub_script_create_file_test.txt.e11979
		  exec_host = node-2-1.local/7
	.
	.
	.
	.
	Job Id: 11980.frontend1.local
		  Job_Name = qsub_script_create_file_test.txt
		  Job_Owner = sripakpa@frontend1.local
		  job_state = Q
		  queue = SP
		  server = frontend1.local
		  Checkpoint = u
		  ctime = Fri Mar 30 21:24:37 2012
		  Error_Path = frontend1.local:/home/sripakpa/QSUB_TEST/qsub_script_create_file_test.txt.e11980
		  Hold_Types = n
	'''

	trial_num=0

	while(True):

		trial_num+=1

		if(trial_num==10): error_exit_with_message("Unsuccessfully requested job_status after 10 trials")

		temp_data_filename="temp_job_stat_V%s.txt" %(str(trial_num).zfill(3) )

		#if(exists(temp_data_filename)): submit_subprocess_allow_retry("rm %s" %(temp_data_filename))

		if(exists(temp_data_filename)): continue

		num_qstat_call=0

		while(True):

			num_qstat_call+=1

			retcode=system( 'qstat -f > %s' %(temp_data_filename))

			sleep(5); #Sleep for a while

			if(num_qstat_call==10): error_exit_with_message("Unsucessfully called qstat for 10 times!")

			if(retcode==0):##Success
				break

		data=open(temp_data_filename, 'r')   #Sept 30, 2010

		break

	#	JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
	job_info_list=[]

	while(True):

		line=data.readline()

		if(line==''): break #End of file!

		if(line[0:7]=="Job Id:"):

			JOB_OWNER='BLAH_BLAH_BLAH'
			JOB_NAME='BLAH_BLAH_BLAH'

			job_info={}
			job_info['JOBID']=line.split()[-1]
			job_info['EXEC_HOST']='UNKNOWN'

			line=data.readline()

			####################################################################################################
			if(line.split()[0]=='Job_Name'): #Parse Job_Name then Job_Owner
				JOB_NAME=line.split()[-1]

				while(True):
					line=data.readline()

					#print "DEBUG line.split()=", line.split()

					if(len(line.split())==0): line=data.readline() #April 22, 2012: Deal with special empty line case (see right below)

					if(len(line.split())!=1): break

					JOB_NAME+=line.split()[0]

					'''
					Job_Name = home_sripakpa_minirosetta_April_22_4X_CHEM_SHIFT_minimize_bench
					mark_run_main_folder_FARFAR_denovo_4X4_R2_Retrotransposon_SLAVE_JOBS_0
		      											<--Empty line here!
		  			Job_Owner = sripakpa@node-9-34.local
					'''

				if(line.split()[0]!='Job_Owner'): error_exit_with_message("line.split()[0]!='Job_Owner' for line=%s" %(line))
				if(line.split()[1]!='='): error_exit_with_message("line.split()[1]!='=' for line=%s" %(line))

				JOB_OWNER=line.split()[-1]


			elif(line.split()[0]=='Job_Owner'):  #Parse Job_Name then Job_Owner

				JOB_OWNER=line.split()[-1]

				line=data.readline()

				if(line.split()[0]!='Job_Name'): error_exit_with_message("line[0:9]!='Job_Owner' for line=%s" %(line))

				JOB_NAME=line.split()[-1]

				while(True):
					line=data.readline()

					#print "DEBUG line.split()=", line.split()

					if(len(line.split())==0): line=data.readline() #April 22, 2012: Deal with special empty line case.

					if(len(line.split())!=1): break

					JOB_NAME+=line.split()[0]

				if(line.split()[1]!='='): error_exit_with_message("line.split()[1]!='=' for line=%s" %(line))

			else:
				error_exit_with_message("Invalid first job status line=%s" %(line))
			####################################################################################################

			job_info['JOB_NAME']=JOB_NAME
			job_info['USER']=JOB_OWNER



			while(True):
				line=data.readline()
				if(line.split()[0]!='job_state'): continue
				job_info['STATE']=line.split()[-1]

				#BSUB RUN, PENDING -> QSUB Q R

				break

			if job_info['STATE'] == 'R':
				while(line):
					line=data.readline()
					if len(line.split())==0:continue
					if (line.split()[0]!='exec_host'): continue
					job_info['EXEC_HOST'] = line.split()[-1].split('.local')[0]
					#BSUB EXEC_HOST --> QSUB exec_host
					break

			job_info_list.append(job_info)

	if(exists(temp_data_filename)):
		system("rm %s" %(temp_data_filename))

	return job_info_list

	#jobs_stat_lines=[]

	#jobs_stat_lines.append("%s   %s   %s   %s " %('JOBID', 'JOB_NAME', 'USER', 'STATE'))

	#for job_info in job_info_list:
		#jobs_stat_lines.append("%s   %s   %s   %s " %(job_info['JOBID'], job_info['JOB_NAME'], job_info['USER'], job_info['STATE']))


	#return jobs_stat_lines

###############################################################################################

def get_num_pending_jobs():
	# Running jobs...
	#lines = popen_and_readlines('bjobs -w | grep " PEND " ', True)
	#lines = popen('qstat | grep " Q " ').readlines()

	job_info_list= get_queued_jobs_status()

	lines=[]
	for job_info in job_info_list:
		if(job_info['STATE']=='Q'): lines.append(job_info['JOB_NAME'])



	jobdirs = glob( "SLAVE_JOBS/*/running.txt" ) #["SLAVE_1", "SLAVE_2"] #choose have been call queued or submitted NOT running.txt
	jobdirs = map( lambda x : dirname( abspath( x) ) , jobdirs )

	jobdirs_active = []
	for line in lines: #["SLAVE/1/", "SLAVE/2'"]

		for jobdir in jobdirs:
			job_tag = jobdir.replace('/','_')[1:]
			if (line==job_tag): #String Comparison
				jobdirs_active.append( jobdir )

	return len( jobdirs_active )

def kick_off_slave_jobs( N_JOBS, nhours = 48 ):

	#return # temporary

	num_pending_jobs = get_num_pending_jobs()
	MAX_PEND = 420
	if ( num_pending_jobs > MAX_PEND ):
		print 'TOO MANY PENDING JOBS ', num_pending_jobs
		return

	print_title_text("Enter kick_off_slave_jobs N_JOBS=%d " %(N_JOBS))

	sys.stdout.flush()
	sys.stderr.flush()

	JOBDIR = 'SLAVE_JOBS'
	if not exists( JOBDIR) :
		command = 'mkdir -p '+JOBDIR
		system( command )
		#system( 'chmod 777 -R '+JOBDIR)

	for n in range( N_JOBS ):
		JOBDIR = 'SLAVE_JOBS/%d' % n
		if not exists( JOBDIR) :
			command = 'mkdir -p '+JOBDIR
			print( command )
			system( command )
			#system( 'chmod 777 -R '+JOBDIR)
		sys.stdout.flush()
		sys.stderr.flush()

	# Jobs may already be running.
	#lines = popen_and_readlines( 'qstat ' , True)
	#lines = popen_and_readlines( 'bjobs -u all -w | grep -v node-4-1.l | grep -v node-9-22.l | grep -v node-8-10.l ' , True)

	job_info_list= get_queued_jobs_status()

	lines=[]
	for job_info in job_info_list:
		lines.append(job_info['JOB_NAME'])

	get_num_pending_jobs

	for n in range( N_JOBS ):

		JOBDIR = 'SLAVE_JOBS/%d' % n

		job_already_running = 0
		for line in lines:
			job_tag = abspath( JOBDIR ).replace('/','_')[1:] #SLAVE/1 -> SLAVE_1
			if line == job_tag: #
				job_already_running = 1
				break #don't submit again! early break!

		if job_already_running:
			print "Already queued: ", JOBDIR
			continue

		running_file = JOBDIR+'/running.txt'
		command = 'echo "running" > '+running_file
		print( command )
		system( command )
		#system( 'chmod 777 '+running_file)

		output = JOBDIR + '/slave_jobs.out'
		system( 'echo "initialized" > '+output )
		#system( 'chmod 777 '+output )

		error = JOBDIR + '/slave_jobs.err'
		system( 'echo "initialized" > '+error )
		#system( 'chmod 777 '+error )

		job_tag = abspath(JOBDIR).replace('/','_')[1:] # /home/frank/SLAVE/1 -> home_frank_SLAVE_1

		#########################

		job_script= SLAVE_EXE + ' ' + JOBDIR

		qsub_submit_file="QSUB_JOB.sh" #"qsub_submit_file"

		QSUB_JOB = open( qsub_submit_file, 'w')

		QSUB_JOB.write( '#!/bin/bash\n\n' )

		QSUB_JOB.write( '#PBS -N %s\n\n' %(job_tag))
		QSUB_JOB.write( '#PBS -o %s\n\n' %(output + "_QSUB" ))
		QSUB_JOB.write( '#PBS -e %s\n\n' %(error + "_QSUB" ))
		QSUB_JOB.write( '#PBS -l walltime=%d:00:00\n\n' % nhours )
		#QSUB_JOB.write( '#PBS -l nodes=1:ppn=1\n\n' )

		QSUB_JOB.write( 'cd $PBS_O_WORKDIR\n\n' ) #This doesn't work if the node submitting the job is not the frontend (i.e. is itself a working_node)

		QSUB_JOB.write( 'cat $PBS_NODEFILE > '+JOBDIR+'/nodefile.txt\n\n' ) # trying to track down problem nodes.

		QSUB_JOB.write( '%s >%s 2>%s\n\n' %(job_script, output, error))

		QSUB_JOB.close()


		############################
		#command = 'bsub -W 140:0 -M 4000000  -R ''hname!=node-4-1'' -o %s -e %s -J %s %s %s ' % (output,error,job_tag,SLAVE_EXE,JOBDIR)
		command = 'qsub -V %s' %(qsub_submit_file)
		print( command )
		system( command )

		rosetta_errfile = '%s/rosetta_errfile_location.txt' % JOBDIR
		system( 'echo "initialized" > '+rosetta_errfile )
		#system( 'chmod 777 '+rosetta_errfile )

		sys.stdout.flush()
		sys.stderr.flush()

	print_title_text("Exit kick_off_slave_jobs N_JOBS=%d " %(N_JOBS))
	return 1

##############################################

def check_is_slave_node_broken(jobdir, slave_errfile, verbose=True):

#	if(verbose): print "enter check_is_slave_node_broken(jobdir=%s)" %(jobdir)

	if(exists(slave_errfile)==False):

		while(exists(slave_errfile)==False):
			print "exists(slave_errfile %s)==False" %(slave_errfile)
			sys.stdout.flush()
			sleep(2)
		sleep(5) # sleep for five more seconds to ensure that slave_errfile is properly closed..

	broken_slave_node=jobdir + "/broken_slave_node.txt"
	if(exists(broken_slave_node)): return True

	#DISABLED, May 10, 2012 -- Rhiju
	if (line_counts(slave_errfile)!=0):

		data =safe_open(slave_errfile, 'r', True)

		for line in data:

			if(line.find( 'createJobTmpDir: Unable to create the job level tmp directory' ) > 0):
				print "check_is_slave_node_broken=True for jobdir=%s: %s " %(jobdir, line)
				submit_subprocess( "echo broken slave node > %s" % ( broken_slave_node ), True )
				return True

	return False


##############################################
def check_output_files( output_files_ ):  #MAKE IT SO THAT IF FOUND ERROR...quit right away

	still_running = 0
	done_running = 0

	for output_file in output_files_:
		still_running += 1
		#print "CHECKING ", output_file

		if not exists( output_file ): continue

		lines = popen_and_readlines("tail -n 20 "+output_file, True)

		for line in lines:
			line=line[:-1]
			if (line == 'Exit job, please check errfile to ensure that there are no errors.' ):
				#print "completed!"
				done_running += 1
				break


	still_running = still_running - done_running
	return still_running

#print "checkaroo: ", still_running, done_running
#            if ( len( line ) > 22 and  line[:23] == "Successfully completed." ):


#IF NO SLAVE IS FREE, this actually waits for the slave to be free before submitting the job
def find_and_submit_job_to_a_slave( qsub_scripts ):

	#print "blah_Blah_blah"
	#sys.stdout.flush()
	#sys.stderr.flush()

	jobdirs = glob( "SLAVE_JOBS/*/running.txt" )
	jobdirs = map( lambda x : dirname( abspath( x) ) , jobdirs )

	num_slave = -1

	# Look at bjobs queue to make sure this job is still running...
	while len( qsub_scripts ) > 0:

		# Running jobs...
		#lines = popen_and_readlines('qstat | grep "R"', True)
		#lines = popen_and_readlines('bjobs -u all -w | grep " RUN " | grep -v node-4-1.l | grep -v node-9-22.l | grep -v node-8-10.l ', True)

		job_info_list= get_queued_jobs_status()

		#print job_info_list
		#sys.stdout.flush()

		lines=[]
		bad_exec_hosts = [ 'node-4-31' ]
		for job_info in job_info_list:
			exec_host = job_info['EXEC_HOST']
			#print exec_host
			if exec_host in bad_exec_hosts: continue
			if(job_info['STATE']=='R'): lines.append(job_info['JOB_NAME'])


#		print lines
#		sys.stdout.flush()
#		sys.stderr.flush()

		jobdirs_active = []
		for line in lines:

 			for jobdir in jobdirs:
				job_tag = jobdir.replace('/','_')[1:]
				if (line==job_tag) > 0: #String Comparison

					# new -- Aug. 2011. Some nodes just aren't running.
					slave_ok_file = jobdir+'/slave_ok.txt'
					if not exists( slave_ok_file ):	continue

					jobdirs_active.append( jobdir )

		#print 'JOBDIRS_ACTIVE: ', jobdirs_active
		#stdout.flush()

		num_slave = -1
		for jobdir in jobdirs_active:  #Check all the slave jobs that are still running/not pending...in case some were killed?

			command_file_name = jobdir + '/run_this_script.txt'

			if not exists( command_file_name ):  #BASICALLY IF COMMAND_FILE_NAME EXIST THEN LAST JOB IS NOT FINISH YET

				####################Check for errors of the previous rosetta job submission#############################
				error_check_method_1=False
				error_check_method_2=False

				slave_errfile= jobdir + '/slave_jobs.err'

				if(check_is_slave_node_broken(jobdir, slave_errfile)): continue

				# Check if node is broken? -- DISABLED, May 10, 2012 -- Rhiju
				if (line_counts(slave_errfile)!=0): #Used to be !=1 before May 4, 2012..
					print "line_counts(slave_errfile)= ", line_counts(slave_errfile)
					error_check_method_1=True

				rosetta_errfile=''
				if(exists( jobdir + "/rosetta_errfile_location.txt" ) and line_counts(jobdir + "/rosetta_errfile_location.txt")!=0):
					rosetta_errfile = popen_and_readlines("tail -n 1 %s" %(jobdir + "/rosetta_errfile_location.txt"), True)[0]
					rosetta_errfile = rosetta_errfile[:-1] #remove the new line (\n) character
					if(exists(rosetta_errfile)):
						if (line_counts(rosetta_errfile)!=0): print 'ERRRRORRRRRRRRRRRRRRRRRRRRRRRR', rosetta_errfile
						#error_check_method_2=True

				if(error_check_method_1 or error_check_method_2):
					print "slave_errfile: " , slave_errfile,
					print "rosetta_errfile: " , rosetta_errfile,
					print "error_check_method_1 ", error_check_method_1
					print "error_check_method_2 ", error_check_method_2
					kill_all_slave_jobs_and_exit()

				####Temporary##############################################################
				qsub_script = qsub_scripts[ 0 ]
				rosetta_errfile_location = qsub_script.replace('.qsub','.err') ###THIS IS HACKY
				submit_subprocess( "echo %s >> %s" % ( rosetta_errfile_location, jobdir + "/rosetta_errfile_location.txt"), True )

				#############################################################################################################
				fid = safe_open( command_file_name, 'w' )
				fid.write( '%s\n' % qsub_script )
				fid.close()

				#system( 'chmod 777 '+command_file_name ) # To allow deletion.

				num_slave = int( basename( jobdir ) ) #basename for '/foo/bar/' returns 'bar', again this relies on the specific form of the job_name..
				assert( exists( command_file_name ) )

				print 'assigning  %s to slave %d' % ( qsub_script, num_slave )
				stdout.flush()

				del( qsub_scripts[ 0 ] )
				if len( qsub_scripts ) == 0: break

		if len( qsub_scripts ) > 0:
			# Did not find anything ... wait a couple seconds
			print "waiting for a slave to free up."
			stdout.flush()
			sleep( 2 )

	return num_slave

def finish_jobs():
	jobdirs = glob( 'SLAVE_JOBS/*/running.txt' )
	jobdirs.sort()
	jobdirs = map( lambda x:dirname( x ), jobdirs )
	for jobdir in jobdirs:
		fid = safe_open( jobdir + '/finished.txt','w' )
		fid.write( ' Finished!! ' )
		fid.close()


def get_current_num_slave_jobs():
	jobdirs = glob( 'SLAVE_JOBS/*/running.txt' )
	return len( jobdirs )

def get_current_num_active_slave_jobs():
	jobdirs = glob( 'SLAVE_JOBS/*/run_this_script.txt' )
	return len( jobdirs )

def kick_off_more_slave_jobs_if_necessary( queue_num, N_JOBS, nhours = 48 ):

	current_num_slave_jobs = get_current_num_slave_jobs()
	if ( current_num_slave_jobs >= N_JOBS ): return

	current_num_active_slave_jobs = get_current_num_active_slave_jobs()
	current_num_free_slave_jobs = current_num_slave_jobs - current_num_active_slave_jobs

	if ( queue_num <= current_num_free_slave_jobs ): return

	needed_jobs = min( N_JOBS, current_num_active_slave_jobs + queue_num )

	print 'Maximum num_jobs: %d.  Current num slave jobs: %d.  Current num active slave jobs: %s. Queue num %d.   Kicking off total of jobs: %d' % \
	      (N_JOBS, current_num_slave_jobs, current_num_active_slave_jobs, queue_num, needed_jobs )

	print 'about to KICK OFF MORE SLAVE JOBS ', needed_jobs
	kick_off_slave_jobs( needed_jobs, nhours )
	sleep( 2 )

##############################################
# Now this script has the power to kick off more slave jobs if it wants to.
#  The limit of the total number of jobs is, of course, N_JOBS.
#

def condor_submit( condor_submit_file_, N_JOBS = 0, nhours = 48 ):
	print_title_text('Run: '+ condor_submit_file_)
	sys.stdout.flush()
	sys.stderr.flush()

	lines = safe_open( condor_submit_file_ ).readlines()
	log = ""
	output = ""
	err = "/dev/null"
	exe = ""
	args = ""
	universe = "vanilla"
	queue_num = 0
	for line in lines:
		if len( line ) > 2:
			cols = string.split( line )
			if cols[0] ==  "executable":
				assert( cols[1] == "=" )
				exe = cols[2]
			elif cols[0] == "arguments":
				assert( cols[1] == "=" )
				args = string.join(cols[2:])
			elif cols[0] == "log":
				assert( cols[1] == "=" )
				log = cols[2]
			elif cols[0] == "output":
				assert( cols[1] == "=" )
				output = cols[2]
			elif cols[0] == "error":
				assert( cols[1] == "=" )
				err = cols[2]
			elif cols[0] == "universe":
				assert( cols[1] == "=" )
				universe = cols[2]
			elif (cols[0]).lower() == "queue":
				if len( cols ) > 1: queue_num = int( cols[1] )

#    if ( exe == "" ) or  ( args == "" ) or (output == "") or (queue_num == 0) or ( log == "" ):
#		print "Error in submitting " + condor_submit_file_
#		assert( False )

	# Stolen from script for PBS.
	output_files_ = []

	kick_off_more_slave_jobs_if_necessary( queue_num, N_JOBS, nhours )
	sys.stdout.flush()

	qsub_scripts = []

	already_done = 1

	for q in range( queue_num ):

		sys.stdout.flush()

		q_string = '%d' % q
		args_new = args.replace('$(Process)',q_string).replace('$(PROCESS)',q_string)
		err_new = err.replace('$(Process)',q_string).replace('$(PROCESS)',q_string)
		output_new = output.replace('$(Process)',q_string).replace('$(PROCESS)',q_string)
		output_files_.append( output_new )

		if exists( output_new ) and ( check_output_files( [output_new] ) == 0 ): # done!
			continue  # this will make sure that we don't kick off any more slave jobs
			#kill_all_slave_jobs_and_exit('Error: The file %s already exist??!!, perhap from a previous job submission' %(output_new) )

		already_done = 0

		# Need to prepare a script for qsub
		qsub_script = err_new.replace('.err','.qsub') ###THIS IS HACKY
		fid = safe_open( qsub_script, 'w' )
		fid.write( '#!/bin/bash\n' )

		fid.write( 'cd %s\n\n' % getcwd() )
		command = "%s %s > %s 2> %s\n\n" % ( exe, args_new, output_new, err_new)
		fid.write( command )

		command = "chmod 777 " +output_new+"\n"
		#fid.write( command )

		command = "chmod 777 " +err_new+"\n"
		#fid.write( command )

		cols = string.split( args_new )
		if cols.count( '-out:file:silent' ):
			outfile = cols[  cols.index('-out:file:silent') +1 ]
			command = "chmod 777 " +outfile+"\n" # to allow for removal and stuff.
			#fid.write( command )

		command = "echo ' ' >> %s \n" % \
			  ( output_new )
		fid.write( command )
		fid.write( command )

		command = "echo 'Exit job, please check errfile to ensure that there are no errors.' >> %s \n" % (output_new )
		fid.write( command )
		fid.close()

		qsub_scripts.append( qsub_script )

        # If this is a scheduler script, just run it.
	#        if universe == "scheduler" and not already_done: ##WHAT THE HELL IS THIS
	#            command = "source "+qsub_script
	#            print( command )
	#            system( command )

	print 'NEED TO RUN: ', qsub_scripts

	num_slave = find_and_submit_job_to_a_slave( qsub_scripts ) #qsub_script is the actual job_command submitted to biox

	#if ( already_done ): output_files_ = []  # send back a blank list to flag that we don't need to post-process! The whole thing is done!

	return output_files_



# Moved this out of SWA_util_2, since it depends on get_queued_jobs_status...
def kill_all_slave_jobs_and_exit(exit_message=""):

	JOBDIR = 'SLAVE_JOBS/'
	job_tag = abspath( JOBDIR ).replace('/','_')[1:]

	#job_tag = home_frank_Knotting_job_SLAVE_JOBS

	job_info_list=get_queued_jobs_status()

	#job_info_list is a array of job_info

	#job_info is the dictionary containing ["JOBID"] and ["Job_name"]

	#Job_name = home_frank_RNA_job_SLAVE_JOBS_1
	#Job_name = home_frank_RNA_job_SLAVE_JOBS_2
	#Job_name = home_frank_RNA_job_SLAVE_JOBS_3
	#Job_name = home_frank_RNA_job_SLAVE_JOBS_4
	#Job_name = home_frank_RNA_job_SLAVE_JOBS_10

	#job_name[:len(job_tag)] converts home_frank_RNA_job_SLAVE_JOBS_1 to home_frank_RNA_job_SLAVE_JOBS

	#Job_name = home_frank_Knottin_job_SLAVE_JOBS_1
	#Job_name = home_frank_Knottin_job_SLAVE_JOBS_2
	#Job_name = home_frank_Knottin_job_SLAVE_JOBS_3
	#Job_name = home_frank_Knottin_job_SLAVE_JOBS_4

	for job_info in job_info_list:

		job_name=job_info["Job_name"]
		job_ID=job_info["JOBID"]

		if(job_name[:len(job_tag)]==job_tag):

			kill_command="qdel %s" %(job_ID)

			print kill_command
			system( kill_command )

			sys.stdout.flush()
			sys.stderr.flush()

	#command = 'qdel %s*' %job_tag
	#print( command )
	#system( command )
#	sleep(10)

	sys.stdout.flush()
	sys.stderr.flush()
	exit_message2="Killed all slave jobs and about to exit master script, " + exit_message
	error_exit_with_message(exit_message2)

