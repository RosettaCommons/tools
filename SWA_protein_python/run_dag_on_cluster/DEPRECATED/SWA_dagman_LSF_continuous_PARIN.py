#!/usr/bin/env python


from sys import argv,exit,stdout
import string
from glob import glob
from os import system,popen,getcwd
from os.path import basename,exists,expanduser,dirname,abspath
import os
import time
from time import sleep
from SWA_util import *

#######################################
# This could actually be a class!
#######################################

###########When is this called#######when SWA_dagman_LSF_continuous is imported?###
PYDIR = get_PYDIR()
assert( exists( PYDIR ) )
SLAVE_EXE = PYDIR + '/SWA_dagman_slave.py'

NUM_CLUSTERER_SLAVE_NODE=6
print "NUM_CLUSTERER_SLAVE_NODE=%d" %(NUM_CLUSTERER_SLAVE_NODE)
####################################################################################
def wait_until_clusterer_slave_nodes_run():

	while(True): #Make sure that at least 3 of the clusterer slave jobs are running!

		bjobs_lines = popen_and_readlines('bjobs -w | grep " RUN " ', True)

		running_clusterer_slave_nodes=[]

		for n in range( NUM_CLUSTERER_SLAVE_NODE ):

			slave_dir = 'SLAVE_JOBS/%d' % n

			slave_tag = abspath( slave_dir ).replace('/','_')

			for line in bjobs_lines:
				if( line.find( slave_tag+' ' ) > 0 ): 
					running_clusterer_slave_nodes.append(slave_tag)
					break

		if(len(running_clusterer_slave_nodes)>=(NUM_CLUSTERER_SLAVE_NODE/2)): #At least half of them are running!
			print "At least half (%d) of the clusterer_slave_nodes are running!" %(NUM_CLUSTERER_SLAVE_NODE/2)
			for n in range(len(running_clusterer_slave_nodes)):
				print "running_clusterer_slave_nodes #%d: %s " %(n+1,  running_clusterer_slave_nodes[n] )
			return
		else:		
			print "Waiting for clusterer_slave_nodes to RUN, SO_FAR %d are running" %(len(running_clusterer_slave_nodes))

		sleep(2)

##########################################################################################


def kick_off_slave_jobs( N_JOBS ):
	print_title_text("Enter kick_off_slave_jobs N_JOBS=%d " %(N_JOBS))

	if(exists('SLAVE_JOBS/')): 	kill_all_slave_jobs_and_exit( 'The folder SLAVE_JOBS/ already exist!') 

	sys.stdout.flush()
	sys.stderr.flush()

	for n in range( N_JOBS ):
		SLAVE_DIR = 'SLAVE_JOBS/%d' % n
		if( exists( SLAVE_DIR)==False ):
			command = 'mkdir -p '+SLAVE_DIR
			print( command )
			submit_subprocess_allow_retry( command, True )
		sys.stdout.flush()
		sys.stderr.flush()

  # ENSURE THAT SLAVE_NODE with same name is NOT ALREADY SUBMITTED!.
	existing_bjobs_lines = popen_and_readlines( 'bjobs -w ' , True)

	for n in range( N_JOBS ):
		
		SLAVE_DIR = 'SLAVE_JOBS/%d' % n

		slave_tag = abspath( SLAVE_DIR ).replace('/','_')

		slave_already_queued = False

		for line in existing_bjobs_lines:
			if(line.find( slave_tag+' ' ) > 0):
				slave_already_queued = True
				break
		
		if (slave_already_queued):	 kill_all_slave_jobs_and_exit("slave_node (%s) is already queued! " %(SLAVE_DIR) )			

		slave_errfile = SLAVE_DIR + '/slave_jobs.err'
		slave_outfile = SLAVE_DIR + '/slave_jobs.out'

		##########################################################################################

		if(n==NUM_CLUSTERER_SLAVE_NODE): wait_until_clusterer_slave_nodes_run()

		if(n<NUM_CLUSTERER_SLAVE_NODE): #HACKY. High memory node for clustering!
			#command = 'bsub -W 140:0 -M 6144000 -R "rusage[mem=6144]" -o %s -e %s -J %s %s -slave_dir %s ' % (slave_outfile, slave_errfile, slave_tag, SLAVE_EXE, SLAVE_DIR)
			command = 'bsub -W 140:0 -M 5000000 -R "rusage[mem=5000]" -o %s -e %s -J %s %s -slave_dir %s ' % (slave_outfile, slave_errfile, slave_tag, SLAVE_EXE, SLAVE_DIR)
		else: #Standard node.
			command = 'bsub -W 140:0 -M 4000000 -o %s -e %s -J %s %s -slave_dir %s ' % (slave_outfile, slave_errfile, slave_tag , SLAVE_EXE, SLAVE_DIR)
			#command = 'bsub -W 140:0 -M 5000000 -R "rusage[mem=5000]" -o %s -e %s -J %s %s -slave_dir %s ' % (slave_outfile, slave_errfile, slave_tag, SLAVE_EXE, SLAVE_DIR)


		print( command )
		submit_subprocess( command, True )
        
		sys.stdout.flush()
		sys.stderr.flush()


	print_title_text("Exit kick_off_slave_jobs N_JOBS=%d " %(N_JOBS))
	return 0 #dummy number, not used for anything.

##############################################

def check_is_slave_node_broken(slave_dir, slave_errfile, slave_outfile, verbose=True):

	######If slave_node was determined to broken is previous call, then just return!#########
	broken_slave_node=slave_dir + "/broken_slave_node.txt"
	if(exists(broken_slave_node)): return True

	#########################################################################################

	if(exists(slave_errfile)==False): ####Wait until slave_errfile exist

		while(exists(slave_errfile)==False): 
			print "exists(slave_errfile %s)==False" %(slave_errfile)
			sys.stdout.flush()
			sleep(2)
		sleep(5) # sleep for five more seconds to ensure that slave_errfile is properly closed..


	if(line_counts(slave_errfile)!=0):

		data =safe_open(slave_errfile, 'r', True)	

		for line in data:	
			if(line.find( 'createJobTmpDir: Unable to create the job level tmp directory' ) > 0):
				print "check_is_slave_node_broken=True for slave_dir=%s: %s " %(slave_dir, line)
				submit_subprocess( "echo broken slave node > %s" % ( broken_slave_node ), True )
				data.close()
				return True
		data.close()


		if(exists(slave_outfile)==False): #Wait until slave_outfile exist

			while(exists(slave_outfile)==False): 
				print "exists(slave_outfile %s)==False" %(slave_outfile)
				sys.stdout.flush()
				sleep(2)
			sleep(5) # sleep for five more seconds to ensure that slave_outfile is properly closed (Aug 10, 2011)..


		if(line_counts(slave_outfile)==0): ###If error_message occur in slave_errfile even before slave_outfile is initialize then BROKEN!##
			print "slave_errfile= %s " %(slave_errfile)
			print "slave_outfile= %s " %(slave_outfile)
			print "check_is_slave_node_broken=True for slave_dir=%s: line_counts(slave_errfile)>0 but line_counts(slave_outfile)==0 " %(slave_dir)
			submit_subprocess( "echo broken slave node > %s" % ( broken_slave_node ), True )
			return True


	##########################################################################################

	return False

##############################################

def slave_job_error_check_method_1(slave_errfile, verbose=False):

	if(verbose): print "enter error_check_method_1(slave_errfile=%s)" %(slave_errfile)

	if(exists(slave_errfile)==False): 
		if(verbose): print "exists(slave_errfile)==False"
		return

	if (line_counts(slave_errfile)!=0):

		exit_message="error_check_method_1=true, line_counts(slave_errfile)= %d, slave_errfile %s " %(line_counts(slave_errfile),  slave_errfile)
		print exit_message

		print "----slave_errfile lines----"
		data =safe_open(slave_errfile, 'r', Is_master=True)	

		line_num=0

		for line in data:	
			line_num+=1
			print "line_num=#", line_num, ": ", line,

		data.close()
		print "----slave_errfile lines----"

		kill_all_slave_jobs_and_exit(exit_message)


##############################################

def error_check_all_slaves(verbose=False):

	slave_dir_list = glob( "SLAVE_JOBS/*/slave_jobs.out" )
	slave_dir_list = map( lambda x : dirname( abspath( x) ) , slave_dir_list )

	for slave_dir in slave_dir_list:  #Check all the slave jobs that are still running/not pending...in case some were killed?

		slave_errfile= slave_dir + '/slave_jobs.err'
		slave_outfile= slave_dir + '/slave_jobs.out'

		if(check_is_slave_node_broken(slave_dir, slave_errfile, slave_outfile)): continue

		slave_job_error_check_method_1(slave_errfile)

	stdout.flush()

##############################################
def check_if_slave_node_disappear(job_info_list):

	bjobs_lines = popen_and_readlines('bjobs -w | grep " RUN " ', True)

	for job in job_info_list:
		found_slave_node=False

		if(job['already_done']==True): continue

		if(job.has_key('slave_node_dir')==False):
			print "ERROR! job: ", job  
			kill_all_slave_jobs_and_exit("job does not have key 'slave_node_dir'!" )	

		slave_tag = job['slave_node_dir'].replace('/','_')

		for line in bjobs_lines:
			if( line.find( slave_tag+' ' ) > 0 ): #OK the +' ' is to distinguish between for example SLAVE_JOBS_19 and SLAVE_JOBS_191
				#print basename(job['slave_node_dir']), " found in line %s " %(line)
				found_slave_node=True
				break

		if(found_slave_node==False): 
			print "-----------------------------------------------------------------------------------------------------------"
			print "ERROR!: slave_node %s is no longer in the bjobs queue!" %(job['slave_node_dir'])
			print "This slave_node was running the job:", job
			print "bjobs -w | grep \" RUN \" before exit:"
			for line in bjobs_lines:
				print line,
			print "-----------------------------------------------------------------------------------------------------------"
			submit_subprocess("bjobs -w > All_bjobs_BEFORE_EXIT.txt")
			kill_all_slave_jobs_and_exit("ERROR! slave_node (%s) is no longer in the bjobs queue!" %(job['slave_node_dir']) )	 


##############################################
def get_run_this_script_file(slave_dir):

	return (slave_dir + "/run_this_script.txt")

##############################################

def ensure_previous_job_has_no_error(slave_dir): ####Check for errors of the previous rosetta job submission

	rosetta_errfile=''

	rosetta_errfile_location=slave_dir + "/rosetta_errfile_location.txt"

	if( ( exists( rosetta_errfile_location ) ) and ( line_counts( rosetta_errfile_location )!=0 ) ):

		rosetta_errfile = popen_and_readlines( "tail -n 1 %s" %( rosetta_errfile_location ), Is_master=True )[0]
		rosetta_errfile = rosetta_errfile[:-1] #remove the new line (\n) character

		if( ( exists(rosetta_errfile) ) and ( line_counts(rosetta_errfile)!=0 ) ): 
			print "rosetta_errfile: " , rosetta_errfile,
			print "error_check_method_2 ", error_check_method_2 
			kill_all_slave_jobs_and_exit()


##############################################

def get_active_slave_dir_list(Is_clusterer_job):

	slave_dir_list = glob( "SLAVE_JOBS/*/slave_jobs.out" )
	slave_dir_list = map( lambda x : dirname( abspath( x ) ) , slave_dir_list )

	slave_dir_list.sort()

	running_bjobs_lines = popen_and_readlines('bjobs -w | grep " RUN " ', True) #FOR checking if slave_node is still running

	active_slave_dir_list = []

	for slave_dir in slave_dir_list:
		
		#####Make sure that slave_node is still RUNNING! (i.e not pending or dead!)####
		slave_tag = slave_dir.replace('/','_')

		slave_is_running=False

		for line in running_bjobs_lines:
			if(line.find( slave_tag + ' ' ) > 0): slave_is_running=True

		if(slave_is_running==False): continue

		#####check that job is not broken THEN make sure that there is no error!#######

		slave_errfile= slave_dir + '/slave_jobs.err'
		slave_outfile= slave_dir + '/slave_jobs.out'

		if(check_is_slave_node_broken(slave_dir, slave_errfile, slave_outfile)): continue

		slave_job_error_check_method_1(slave_errfile, verbose=False)

		####Aug 12, 2011####Check that slave_jobs.out is not empty.##########

		if(exists(slave_outfile)==False):  kill_all_slave_jobs_and_exit("slave_outfile (%s) doesn't exist!!" %(slave_outfile) )

		if(line_counts(slave_outfile)==0): 
			print "Skipping slave_outfile (%s) since line_counts(slave_outfile)==0)" %(slave_outfile) 
			continue

		#####See if slave_node is still running a prior job#############################################
	
		run_this_script_file = get_run_this_script_file(slave_dir)

		if (exists( run_this_script_file ) ): continue #BASICALLY IF run_this_script_file EXIST THEN LAST JOB IS NOT FINISH YET

		######HACKY, Allocate special high RAM memory slave_node to clusterer###########################

		try:
			slave_ID=int(slave_dir.split('/')[-1])
		except:
			kill_all_slave_jobs_and_exit("Cannot get slave_ID from slave_dir (%s)" %(slave_dir))

		if( (Is_clusterer_job==False) and (slave_ID<NUM_CLUSTERER_SLAVE_NODE ) ): continue

		if( (Is_clusterer_job==True ) and (slave_ID>=NUM_CLUSTERER_SLAVE_NODE) ): continue

		######################################################################

		ensure_previous_job_has_no_error(slave_dir)

		active_slave_dir_list.append( slave_dir )

	return active_slave_dir_list

##############################################
def get_next_job_queue_ID(job_queue_ID, job_info_list):

	if(job_queue_ID>=len(job_info_list)): kill_all_slave_jobs_and_exit("job_queue_ID=(%s)>=(%s)=len(job_info_list)" %(job_queue_ID,len(job_info_list) ) )

	while(True): #find next job in the list that is not yet done.

		job_queue_ID+=1

		if(job_queue_ID==len(job_info_list)): break #break if reach end of list.

		if(job_info_list[job_queue_ID]['already_done']==False): break		 #break once find a job that not yet done.			
		
	return job_queue_ID


##############################################

def submit_job(job_info, slave_dir):

	job_script_filename =job_info['job_script_filename'] 
	err_log_filename    =job_info['err_log_filename']
	done_signal_filename=job_info['done_signal_filename']	

	submit_subprocess( "echo %s >> %s" % ( err_log_filename, slave_dir + "/rosetta_errfile_location.txt"), True )

	run_this_script_file = get_run_this_script_file(slave_dir)

	fid = safe_open( run_this_script_file, 'w' )
	fid.write( '%s\n' % job_script_filename )
	fid.write( '%s\n' % done_signal_filename )
	fid.close()

	if( exists( run_this_script_file ) == False ): kill_all_slave_jobs_and_exit("run_this_script_file (%s) doesn't exist!" %(run_this_script_file) )

	print 'assigning %s to slave %d' % ( job_script_filename, int( basename( slave_dir ) ) )

	stdout.flush()

##############################################

def submit_jobs_to_slave( job_info_list, Is_clusterer_job):

	verbose=True

	num_act_jobs_submitted=0

	job_queue_ID=-1

	while( True ):

		active_slave_dir_list=get_active_slave_dir_list(Is_clusterer_job) #Get active jobs that are not pending, not dead and not broken!

		if(len(active_slave_dir_list)==0):
			print "waiting for slave nodes to free up."
			sleep( 2 )

		for slave_dir in active_slave_dir_list:

			job_queue_ID=get_next_job_queue_ID(job_queue_ID, job_info_list) #Get the next job_queue_ID that is not yet done!

			############################################
			if(job_queue_ID==len(job_info_list)): #FINISH SUBMITTING ALL JOBS!!!!

				if(num_act_jobs_submitted==0):
					print "WARNING: job_queue_ID==len(job_info_list), even before submitting any job!, len(job_info_list)=%s" %( len(job_info_list) )

				stdout.flush()
				return job_info_list
			############################################

			submit_job(job_info_list[job_queue_ID], slave_dir)

			job_info_list[job_queue_ID]['slave_node_dir']=slave_dir
			
			num_act_jobs_submitted+=1

			
##############################################

def send_finish_signal_to_slave_nodes():

	slave_dir_list = glob( 'SLAVE_JOBS/*/slave_jobs.out' )

	slave_dir_list.sort()

	slave_dir_list = map( lambda x:dirname( x ), slave_dir_list )

	for slave_dir in slave_dir_list:

		fid = safe_open( slave_dir + '/finished.txt','w' )
		fid.write( ' Finished!! ' )
		fid.close()


##############################################
def delete_log_files(job_info_list, keep_some_log_file):

	
	for n in range(len(job_info_list)):

		job_info=job_info_list[n]

		out_log_filename=job_info['out_log_filename']
		err_log_filename=job_info['err_log_filename']
		job_script_filename=job_info['job_script_filename']
		done_signal_filename=job_info['done_signal_filename']

		if(keep_some_log_file==True): 
			if(out_log_filename.count('KEEP_LOG_FILE/')>0): continue

		command='rm %s ' %(out_log_filename)
		if(exists(out_log_filename)): 
			print command
			submit_subprocess_allow_retry(command,  True )

		command='rm %s ' %(err_log_filename)
		if(exists(err_log_filename)): 
			print command
			submit_subprocess_allow_retry(command,  True )

		command='rm %s ' %(job_script_filename)
		if(exists(job_script_filename)): 
			print command
			submit_subprocess_allow_retry(command,  True )

		command='rm %s ' %(done_signal_filename)
		if(exists(done_signal_filename)): 
			print command
			submit_subprocess_allow_retry(command,  True )

#		command='rm %s ' %(done_signal_filename+"'\n'")
#		if(exists(done_signal_filename+'\n')): 
#			print command
#			submit_subprocess_allow_retry(command,  True )


##############################################
def condor_submit( condor_submit_file_ , reducer_command):

	condor_submit_start_time=time.time()


	print_title_text('Run: '+ condor_submit_file_)
	sys.stdout.flush()
	sys.stderr.flush()

	lines = safe_open( condor_submit_file_ ).readlines()


	exe= ""
	args= ""
	log_foldername= ""
	mapper_infiles= ""
	mapper_outfiles= ""
	reducer_infiles= ""
	reducer_outfiles= ""
#	reducer_rm_files= ""
	for line in lines:
		if len( line ) > 2: 
			cols = string.split( line )
			if cols[0] ==  "executable":
				assert( cols[1] == "=" )
				exe = cols[2]
			elif cols[0] == "arguments":
				assert( cols[1] == "=" )
				args = string.join(cols[2:])
			elif cols[0] == "log_foldername":
				assert( cols[1] == "=" )
				log_foldername = cols[2]
			elif cols[0] == "mapper_infiles": #USED ONLY IN SWA2
				assert( cols[1] == "=" )
				mapper_infiles = string.join(cols[2:])
			elif cols[0] == "mapper_outfiles":
				assert( cols[1] == "=" )
				mapper_outfiles = string.join(cols[2:])
			elif cols[0] == "reducer_infiles": #USED ONLY IN SWA2
				assert( cols[1] == "=" )
				reducer_infiles = string.join(cols[2:])
			elif cols[0] == "reducer_outfiles": 
				assert( cols[1] == "=" )
				reducer_outfiles = string.join(cols[2:])
			elif (cols[0]).lower() == "queue":
				if len( cols ) > 1: queue_num = int( cols[1] )

	reducer_outfile_list=reducer_outfiles.split() 
	for reducer_outfile in reducer_outfile_list:
		if(exists(reducer_outfile)): kill_all_slave_jobs_and_exit("reducer_outfile %s already exist before job submission!" %(reducer_outfile) )	


	# Stolen from script for PBS.
	all_job_info_list = []
	#actually_submit_job_info_list = []

	for q in range( queue_num + 1 ):

		reducer_job=False

		if(q==queue_num): #reducer_job
			reducer_job=True
			job_name='reducer'
			if(reducer_command==""): continue
		else:
			reducer_job=False
			job_name='%d' % q #mapper_job


		job_info = get_job_info(log_foldername, job_name)

		out_log_filename=job_info['out_log_filename']
		err_log_filename=job_info['err_log_filename']
		job_script_filename=job_info['job_script_filename']
		done_signal_filename=job_info['done_signal_filename']

		outfiles=""
		if(reducer_job):
			outfiles=reducer_outfiles
		else:
			outfiles=mapper_outfiles.replace('$(Process)',job_name)

		outfile_list=outfiles.split() 

		########CHECK IF JOB HAVE PREVIOUSLY BEEN SUCCESSFULLY COMPLETED###################
		if exists( done_signal_filename ):

			already_done=True

			if(reducer_job): kill_all_slave_jobs_and_exit("done_signal for reducer_job %s already exist before job submission!" %(done_signal_filename) )
		
			for outfile in outfile_list:
				if(exists(outfile)==False):	print "WARNING: outfile %s doesn't exist" %(outfile)				
#				if(exists(outfile)==False):	kill_all_slave_jobs_and_exit("outfile %s doesn't exist. ALL outfiles= %s" %(outfile, outfiles) )				

				print "job %s is already successfully_completed " %(outfile)
			

		else:
			#######MAKE SURE THAT outfiles and logs files from previously failed runs are deleted####
			already_done=False

			if(exists(out_log_filename)):
				print ("job %s WAS NOT successfully_completed. Removing out_log_filename = %s " %(out_log_filename, out_log_filename)) 
				submit_subprocess( 'rm %s' %(out_log_filename), True)

			if(exists(err_log_filename)):
				print ("job %s WAS NOT successfully_completed. Removing err_log_filename = %s " %(out_log_filename, err_log_filename)) 
				submit_subprocess( 'rm %s' %(err_log_filename), True)
		
			if(exists(job_script_filename)):
				print ("job %s WAS NOT successfully_completed. Removing job_script_filename = %s " %(out_log_filename, job_script_filename)) 
				submit_subprocess( 'rm %s' %(job_script_filename) , True)
		
			for outfile in outfile_list:

				if(exists(outfile)): 
					print ("job %s WAS NOT successfully_completed.	 Removing outfile = %s " %(out_log_filename, outfile)) 
					submit_subprocess( 'rm %s' %(outfile), True )

			######CREATE THE JOB_SCRIPT FOR THIS JOB #######################################################


			local_safe_mkdir( log_foldername , True)
			local_safe_mkdir( dirname(out_log_filename) , True)
			local_safe_mkdir( dirname(err_log_filename) , True)
			local_safe_mkdir( dirname(job_script_filename) , True)
			local_safe_mkdir( dirname(done_signal_filename) , True)

			for outfile in outfile_list:
				local_safe_mkdir( dirname(outfile) , True)

			fid = safe_open( job_script_filename, 'w' )

			if(reducer_job):
				fid.write('wrapper_script = SWA_reducer.py\n')
				fid.write('command = %s\n' % reducer_command )
				fid.write('job_name = %s\n' % (job_name) )
				fid.write('log_foldername = %s\n' % (log_foldername) )
				fid.write('num_mapper_jobs %d\n' % (queue_num) )

			else:
				args_actual = args.replace('$(Process)',job_name)

				fid.write( '#!/bin/bash\n' )

				fid.write( 'cd %s\n\n' % getcwd() )
				command = "%s %s > %s 2> %s\n\n" % ( exe, args_actual, out_log_filename, err_log_filename) 
				fid.write( command )

			fid.close()

		job_info['already_done']=already_done

		all_job_info_list.append( job_info )


	Is_clusterer_job=False

	if(condor_submit_file_.count('cluster')>0): 
		Is_clusterer_job=True
		print "%s is a clusterer job!" %(condor_submit_file_)
	else:
		Is_clusterer_job=False
		print "%s is NOT a clusterer job!" %(condor_submit_file_)

	all_job_info_list=submit_jobs_to_slave( all_job_info_list, Is_clusterer_job )



#	if(queue_num!=0):		
#	else:
#		if(len(all_job_info_list)!=0): kill_all_slave_jobs_and_exit("queue_num==0 but len(all_job_info_list)!=0")
#		#consistency check
#		print "Queue_num=0 for JOB %s, no actual job submission!" %(condor_submit_file_)


	print "A total of %f seconds is used to submit the JOB: %s " %( time.time()- condor_submit_start_time, condor_submit_file_)

	return all_job_info_list


