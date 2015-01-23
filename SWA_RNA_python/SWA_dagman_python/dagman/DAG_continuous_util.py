#!/usr/bin/env python


######################################################################
from DAG_general_util import get_job_info, check_if_job_in_done


######################################################################
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.utility.master_util import *


from SWA_dagman_python.scheduler.scheduler_util import *

######################################################################

from sys import argv,exit,stdout
import string
from glob import glob
from os import system,popen,getcwd
from os.path import basename,expanduser,dirname,abspath
import os
import time
from time import sleep


######################################################################


#######################################
# This could actually be a class!
#######################################

SEPERATE_CLUSTERER_SLAVE_NODES=False
NUM_CLUSTERER_SLAVE_NODE=-1 #8
print "SEPERATE_CLUSTERER_SLAVE_NODES=%s | NUM_CLUSTERER_SLAVE_NODE=%d" %(SEPERATE_CLUSTERER_SLAVE_NODES, NUM_CLUSTERER_SLAVE_NODE)
if(SEPERATE_CLUSTERER_SLAVE_NODES==True and NUM_CLUSTERER_SLAVE_NODE<=0): error_exit_with_message("SEPERATE_CLUSTERER_SLAVE_NODES==True and NUM_CLUSTERER_SLAVE_NODE<=0")


####################################################################################
def wait_until_clusterer_slave_nodes_run():

	if(SEPERATE_CLUSTERER_SLAVE_NODES==False): error_exit_with_message("SEPERATE_CLUSTERER_SLAVE_NODES==False!")

	while(True): #Make sure that at least 3 of the clusterer slave jobs are running!

		running_job_tag_list=get_queued_job_name_list("RUN", ignore_problem_nodes=True)

		running_clusterer_slave_nodes=[]

		for n in range( NUM_CLUSTERER_SLAVE_NODE ):

			slave_dir = 'SLAVE_JOBS/%d' % n

			slave_tag = abspath( slave_dir ).replace('/','_')[1:]

			if(running_job_tag_list.count(slave_tag)>0): running_clusterer_slave_nodes.append(slave_tag)

		if(len(running_clusterer_slave_nodes)>=(NUM_CLUSTERER_SLAVE_NODE/2)): #At least half of them are running!
			print "At least half (%d) of the clusterer_slave_nodes are running!" %(NUM_CLUSTERER_SLAVE_NODE/2)
			for n in range(len(running_clusterer_slave_nodes)):
				print "running_clusterer_slave_nodes #%d: %s " %(n+1,  running_clusterer_slave_nodes[n] )
			return
		else:
			print "Waiting for clusterer_slave_nodes to RUN, SO_FAR %d are running" %(len(running_clusterer_slave_nodes))

		sleep(2)

##########################################################################################
def assert_slave_jobs_not_already_queued(N_JOBS):

	print_title_text("Enter assert_slave_jobs_not_already_queued | N_JOBS=%d " %(N_JOBS))

  # ENSURE THAT SLAVE_NODE with same name is NOT ALREADY SUBMITTED!.
	existing_job_tag_list= get_queued_job_name_list("ALL", ignore_problem_nodes=True)

	for job_num in range( N_JOBS ):

		print "job_num=%s" %(job_num)

		SLAVE_DIR = 'SLAVE_JOBS/%d' %(job_num)

		slave_tag = abspath( SLAVE_DIR ).replace('/','_')[1:]

		slave_already_queued = False

		if(existing_job_tag_list.count(slave_tag)>0): slave_already_queued = True

		if(slave_already_queued):	error_exit_with_message("slave_node (%s) is already queued! " %(SLAVE_DIR) )

	print_title_text("Exit assert_slave_jobs_not_already_queued | N_JOBS=%d " %(N_JOBS))

##########################################################################################

def MPI_DAG_slave_wrapper(slave_num, TOTAL_SLAVE_JOBS, SLAVE_EXE): #This is the only function in this .py file that is called by the slave and not the master!

	print_title_text("Enter MPI_DAG_slave_wrapper | slave_num=%d of TOTAL_SLAVE_JOBS=%d " %(slave_num, TOTAL_SLAVE_JOBS))

	SLAVE_DIR = 'SLAVE_JOBS/%d' %(slave_num)

	if( exists( SLAVE_DIR) ): error_exit_with_message("SLAVE_DIR (%s) already exist!" %(SLAVE_DIR))

	command = 'mkdir -p %s' %(SLAVE_DIR)
	print( command )
	submit_subprocess_allow_retry( command )

	##############################################################
	sys.stdout.flush()
	sys.stderr.flush()

	slave_tag = abspath( SLAVE_DIR ).replace('/','_')[1:]

	slave_outfile = SLAVE_DIR + '/slave_jobs.out' #Commented out on March 31, 2012 #Switch back to this on April 21, 2012
	slave_errfile = SLAVE_DIR + '/slave_jobs.err' #Commented out on March 31, 2012 #Switch back to this on April 21, 2012

	if(exists(slave_outfile)): submit_subprocess("rm %s" %(slave_outfile))
	if(exists(slave_errfile)): submit_subprocess("rm %s" %(slave_errfile))

	outfile= safe_open(slave_outfile, "w")
	errfile= safe_open(slave_errfile, "w")

	os.dup2(outfile.fileno(), sys.stdout.fileno())
	os.dup2(errfile.fileno(), sys.stderr.fileno())

	sys.stdout.flush()
	sys.stderr.flush()
	##############################################################

	slave_job_script= "%s -slave_dir %s " %(SLAVE_EXE, SLAVE_DIR)

	if(SEPERATE_CLUSTERER_SLAVE_NODES): error_exit_with_message("SEPERATE_CLUSTERER_SLAVE_NODES not yet implemented in MPI_mode!")

	#Request of wall_time and memory_limit are already done at the MPI qsub submission phase.

	print "slave_job_script=%s" %(slave_job_script)

	submit_subprocess(slave_job_script)

	sys.stdout.flush()
	sys.stderr.flush()

	###########################################################################################################################################################

	print_title_text("Exit MPI_DAG_slave_wrapper | slave_num=%d of TOTAL_SLAVE_JOBS=%d " %(slave_num, TOTAL_SLAVE_JOBS))
	print

##########################################################################################
def kick_off_slave_jobs( N_JOBS, SLAVE_EXE, wall_time ):
	print_title_text("Enter kick_off_slave_jobs | N_JOBS=%d " %(N_JOBS))

	if(exists('SLAVE_JOBS/')): error_exit_with_message( 'The folder SLAVE_JOBS/ already exist!')

	sys.stdout.flush()
	sys.stderr.flush()

	for job_num in range( N_JOBS ):
		SLAVE_DIR = 'SLAVE_JOBS/%d' %(job_num)
		if( exists( SLAVE_DIR)==False ):
			command = 'mkdir -p %s' %(SLAVE_DIR)
			print( command )
			master_submit_subprocess_allow_retry( command )
		sys.stdout.flush()
		sys.stderr.flush()


	for job_num in range( N_JOBS ):

		SLAVE_DIR = 'SLAVE_JOBS/%d' %(job_num)

		slave_tag = abspath( SLAVE_DIR ).replace('/','_')[1:]

		slave_errfile = SLAVE_DIR + '/slave_jobs.err' #Commented out on March 31, 2012 #Switch back to this on April 21, 2012
		slave_outfile = SLAVE_DIR + '/slave_jobs.out' #Commented out on March 31, 2012 #Switch back to this on April 21, 2012

		#slave_errfile = 'slave_jobs.err' #March 31, 2012 #Commented out on April 21, 2012
		#slave_outfile = 'slave_jobs.out' #March 31, 2012 #Commented out on April 21, 2012

		slave_job_script= "%s -slave_dir %s " %(SLAVE_EXE, SLAVE_DIR)
		##########################################################################################

		if(SEPERATE_CLUSTERER_SLAVE_NODES and (job_num==NUM_CLUSTERER_SLAVE_NODE) ): wait_until_clusterer_slave_nodes_run()

		if(SEPERATE_CLUSTERER_SLAVE_NODES and (job_num<NUM_CLUSTERER_SLAVE_NODE ) ): #HACKY. High memory node for clustering!
			queue_job(slave_tag, slave_outfile, slave_errfile, slave_job_script, SLAVE_DIR, walltime=wall_time, memory_limit=8192, memory_reserve=8192)
			#Before April 27, 2012: walltime=140
		else: #Standard node.
			#queue_job(slave_tag, slave_outfile, slave_errfile, slave_job_script, SLAVE_DIR, walltime=wall_time, memory_limit=4096, memory_reserve=0, verbose=True)

			queue_job(slave_tag, slave_outfile, slave_errfile, slave_job_script, SLAVE_DIR, walltime=wall_time, memory_limit=2048, memory_reserve=2048) #May 11, 2012


		sys.stdout.flush()
		sys.stderr.flush()

	##############################Now identify slave nodes by the JOB_NAME. Could switch to identify by the JOB_ID############################################

	all_slave_job_IDs_filename="ALL_SLAVE_JOB_IDS.txt"

	if(exists(all_slave_job_IDs_filename)): submit_subprocess("rm %s" %(all_slave_job_IDs_filename))

	ALL_SLAVE_JOB_IDS = open( all_slave_job_IDs_filename, 'w')
	print "CURRENT_DIRECTORY=%s" %(os.path.abspath("."))
	sleep(5)
	file_count=0

	for job_num in range( N_JOBS ):
		slave_job_ID_filename = "SLAVE_JOBS/%d/JOB_ID.txt" %(job_num)

		if(exists(slave_job_ID_filename)==False): continue #For the case where the slave node is broken.
		file_count+=1

		slave_job_ID_lines = master_safe_open(slave_job_ID_filename).readlines()

		if(len(slave_job_ID_lines)!=1): error_exit_with_message("len(slave_job_ID_lines)=(%s)!=1 for slave_job_ID_filename=(%s)" %(len(slave_job_ID_lines), slave_job_ID_filename))

		ALL_SLAVE_JOB_IDS.write('%s' %(slave_job_ID_lines[0]))

	ALL_SLAVE_JOB_IDS.close()
	print "%s out of %s slave_job_ID_filename found!" %(file_count, N_JOBS)

	###########################################################################################################################################################

	print_title_text("Exit kick_off_slave_jobs | N_JOBS=%d " %(N_JOBS))
	print

	return 0 #dummy number, not used for anything.



##############################################

def send_finish_signal_to_slave_nodes(N_JOBS):

	print_title_text("send_finish_signal_to_slave_nodes N_JOBS=%d " %(N_JOBS))

	if(exists('SLAVE_JOBS/')==False):	master_kill_all_slave_jobs_and_exit( "The folder SLAVE_JOBS/ doesn't exist!")

	sys.stdout.flush()
	sys.stderr.flush()

	for n in range( N_JOBS ):

		SLAVE_DIR = 'SLAVE_JOBS/%d' % n

		if(exists(SLAVE_DIR)==False): master_kill_all_slave_jobs_and_exit( "The folder SLAVE_DIR (%s) doesn't exist!" %(SLAVE_DIR) )

		fid = master_safe_open( SLAVE_DIR + '/finished.txt','w' )
		fid.write( ' Finished!! ' )
		fid.close()

	'''##########OLD WAY COMMENTED OUT ON DEC 14, 2011##########
	slave_dir_list = glob( 'SLAVE_JOBS/*/slave_jobs.out' )

	slave_dir_list.sort()

	slave_dir_list = map( lambda x:dirname( x ), slave_dir_list )

	for slave_dir in slave_dir_list:

		fid = master_safe_open( slave_dir + '/finished.txt','w' )
		fid.write( ' Finished!! ' )
		fid.close()
	'''

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


	if(master_line_counts(slave_errfile)!=0):

		data =master_safe_open(slave_errfile, 'r')

		for line in data:
			if(line.find( 'createJobTmpDir: Unable to create the job level tmp directory' ) > 0):
				print "check_is_slave_node_broken=True for slave_dir=%s: %s " %(slave_dir, line)
				master_submit_subprocess( "echo broken slave node > %s" % ( broken_slave_node ) )
				data.close()
				return True
		data.close()


		if(exists(slave_outfile)==False): #Wait until slave_outfile exist

			while(exists(slave_outfile)==False):
				print "exists(slave_outfile %s)==False" %(slave_outfile)
				sys.stdout.flush()
				sleep(2)
			sleep(5) # sleep for five more seconds to ensure that slave_outfile is properly closed (Aug 10, 2011)..


		if(master_line_counts(slave_outfile)==0): ###If error_message occur in slave_errfile even before slave_outfile is initialize then BROKEN!##
			print "slave_errfile= %s " %(slave_errfile)
			print "slave_outfile= %s " %(slave_outfile)
			print "check_is_slave_node_broken=True for slave_dir=%s: master_line_counts(slave_errfile)>0 but master_line_counts(slave_outfile)==0 " %(slave_dir)
			master_submit_subprocess( "echo broken slave node > %s" % ( broken_slave_node ) )
			return True


	##########################################################################################

	return False

##############################################

def slave_job_error_check_method_1(slave_errfile, verbose=False):

	if(verbose): print "enter error_check_method_1(slave_errfile=%s)" %(slave_errfile)

	if(exists(slave_errfile)==False):
		if(verbose): print "exists(slave_errfile)==False"
		return

	if(master_line_counts(slave_errfile)!=0):

		exit_message="error_check_method_1=true, master_line_counts(slave_errfile)= %d, slave_errfile %s " %(master_line_counts(slave_errfile),  slave_errfile)
		print exit_message

		print "----slave_errfile lines----"
		data =master_safe_open(slave_errfile, 'r')

		line_num=0

		for line in data:
			line_num+=1
			print "line_num=#", line_num, ": ", line,

		data.close()
		print "----slave_errfile lines----"

		master_kill_all_slave_jobs_and_exit(exit_message)


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

	running_job_tag_list = get_queued_job_name_list("RUN", ignore_problem_nodes=True)

	for job in job_info_list:
		found_slave_node=False

		if(job['already_done']==True): continue

		if(job.has_key('slave_node_dir')==False):
			print "ERROR! job: ", job
			master_kill_all_slave_jobs_and_exit("job does not have key 'slave_node_dir'!" )

		slave_tag = job['slave_node_dir'].replace('/','_')[1:]

		if(running_job_tag_list.count(slave_tag)>0):
			#print basename(job['slave_node_dir']), " found in line %s " %(line)
			found_slave_node=True


		if(found_slave_node==False):
			print "-----------------------------------------------------------------------------------------------------------"
			print "ERROR!: slave_tag (%s) is no longer in the running_job_tag_list!" %(job['slave_node_dir'])
			print "This slave_node was running the job:", job
			print "queued_jobs_stat_lines BEFORE EXIT:"

			QUEUED_JOBS_OUTFILE = open( "All_QUEUED_JOBS_BEFORE_EXIT.txt", 'w')
			queued_jobs_stat_lines=get_queued_jobs_status_lines()

			for n in range(len(queued_jobs_stat_lines)):
				print "#%4d: %s" %(n, queued_jobs_stat_lines[n])
				QUEUED_JOBS_OUTFILE.write('#%4d: %s\n' %(n, queued_jobs_stat_lines[n]))
			print "-----------------------------------------------------------------------------------------------------------"
			QUEUED_JOBS_OUTFILE.close()

			master_kill_all_slave_jobs_and_exit("ERROR! slave_node (%s) is no longer in the running_job_tag_list!" %(job['slave_node_dir']) )


##############################################
def get_run_this_script_file(slave_dir):

	return (slave_dir + "/run_this_script.txt")

##############################################

def ensure_previous_job_has_no_error(slave_dir): ####Check for errors of the previous rosetta job submission

	rosetta_errfile=''

	rosetta_errfile_location=slave_dir + "/rosetta_errfile_location.txt"

	if( ( exists( rosetta_errfile_location ) ) and ( master_line_counts( rosetta_errfile_location )!=0 ) ):

		rosetta_errfile = master_popen_and_readlines( "tail -n 1 %s" %( rosetta_errfile_location ) )[0]
		rosetta_errfile = rosetta_errfile[:-1] #remove the new line (\n) character

		if( ( exists(rosetta_errfile) ) and ( master_line_counts(rosetta_errfile)!=0 ) ):
			if is_valid_error(rosetta_errfile):
				print "NON-EMPTY rosetta_errfile (%s)!" %(rosetta_errfile)
				print "Error slave_node=%s " %(slave_dir)
				master_kill_all_slave_jobs_and_exit()

##############################################

def is_valid_error(errfile):
	for line in errfile:
	####if '<TEXT THAT SHOULD NOT BE CAUGHT AS AN ERROR>' in line: continue
		if 'HESSIN for (i,i):' in line:	continue
		if 'G for (i):' in line:		continue
		if not len(line):				continue
		return True
	return False

##############################################

def get_active_slave_dir_list(Is_clusterer_job):

	slave_dir_list = glob( "SLAVE_JOBS/*/slave_jobs.out" )
	slave_dir_list = map( lambda x : dirname( abspath( x ) ) , slave_dir_list )

	slave_dir_list.sort()

	running_job_tag_list = get_queued_job_name_list("RUN", ignore_problem_nodes=True)

	#May 09, 2012: Hacky for debugging purpose!
	update_scheduler_queue_log()

	active_slave_dir_list = []

	for slave_dir in slave_dir_list:

		#####Make sure that slave_node is still RUNNING! (i.e not pending or dead!)####
		slave_tag = slave_dir.replace('/','_')[1:]

		slave_is_running=False

		if(running_job_tag_list.count(slave_tag)>0):
			slave_is_running=True

		if(slave_is_running==False): continue

		#####check that job is not broken THEN make sure that there is no error!#######

		slave_errfile= slave_dir + '/slave_jobs.err'
		slave_outfile= slave_dir + '/slave_jobs.out'

		if(check_is_slave_node_broken(slave_dir, slave_errfile, slave_outfile)): continue

		slave_job_error_check_method_1(slave_errfile, verbose=False)

		####Aug 12, 2011####Check that slave_jobs.out is not empty.##########

		if(exists(slave_outfile)==False):  master_kill_all_slave_jobs_and_exit("slave_outfile (%s) doesn't exist!!" %(slave_outfile) )

		if(master_line_counts(slave_outfile)==0):
			print "Skipping slave_outfile (%s) since master_line_counts(slave_outfile)==0)" %(slave_outfile)
			continue

		#####See if slave_node is still running a prior job#############################################

		run_this_script_file = get_run_this_script_file(slave_dir)

		if (exists( run_this_script_file ) ): continue #BASICALLY IF run_this_script_file EXIST THEN LAST JOB IS NOT FINISH YET

		######HACKY, Allocate special high RAM memory slave_node to clusterer###########################

		try:
			slave_ID=int(slave_dir.split('/')[-1])
		except:
			master_kill_all_slave_jobs_and_exit("Cannot get slave_ID from slave_dir (%s)" %(slave_dir))

		if(SEPERATE_CLUSTERER_SLAVE_NODES):
			if( (Is_clusterer_job==False) and (slave_ID<NUM_CLUSTERER_SLAVE_NODE ) ): continue
			if( (Is_clusterer_job==True ) and (slave_ID>=NUM_CLUSTERER_SLAVE_NODE) ): continue

		######################################################################

		ensure_previous_job_has_no_error(slave_dir)

		active_slave_dir_list.append( slave_dir )

	return active_slave_dir_list

##############################################
def get_next_job_queue_ID(job_queue_ID, job_info_list):

	if(job_queue_ID>=len(job_info_list)): master_kill_all_slave_jobs_and_exit("job_queue_ID=(%s)>=(%s)=len(job_info_list)" %(job_queue_ID,len(job_info_list) ) )

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

	master_submit_subprocess( "echo %s >> %s" % ( err_log_filename, slave_dir + "/rosetta_errfile_location.txt") )

	run_this_script_file = get_run_this_script_file(slave_dir)

	fid = master_safe_open( run_this_script_file, 'w' )
	fid.write( '%s\n' % job_script_filename )
	fid.write( '%s\n' % done_signal_filename )
	fid.close()

	if( exists( run_this_script_file ) == False ): master_kill_all_slave_jobs_and_exit("run_this_script_file (%s) doesn't exist!" %(run_this_script_file) )

	print 'assigning %s to slave %d' % ( job_script_filename, int( basename( slave_dir ) ) )

	stdout.flush()

##############################################

def submit_jobs_to_slave( job_info_list, Is_clusterer_job):

	verbose=True

	num_act_jobs_submitted=0

	job_queue_ID=-1

	#hack_count=0

	while( True ):

		error_check_all_slaves() #Add on April 21, 2012

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

		#hack_count+=1

		#if(hack_count==10):
		#	print "hack_count==10!, master_kill_all_slave_jobs_and_exit()"
		#	master_kill_all_slave_jobs_and_exit("hack_count==10")

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
			master_submit_subprocess_allow_retry(command)

		command='rm %s ' %(err_log_filename)
		if(exists(err_log_filename)):
			print command
			master_submit_subprocess_allow_retry(command)

		command='rm %s ' %(job_script_filename)
		if(exists(job_script_filename)):
			print command
			master_submit_subprocess_allow_retry(command)

		command='rm %s ' %(done_signal_filename)
		if(exists(done_signal_filename)):
			print command
			master_submit_subprocess_allow_retry(command)

#		command='rm %s ' %(done_signal_filename+"'\n'")
#		if(exists(done_signal_filename+'\n')):
#			print command
#			master_submit_subprocess_allow_retry(command)


##############################################
def DAG_submit( DAG_submit_file_ , reducer_command):

	DAG_submit_start_time=time.time()

	print_title_text('Run: '+ DAG_submit_file_)
	sys.stdout.flush()
	sys.stderr.flush()


	lines = master_safe_open( DAG_submit_file_ ).readlines()
	lines = [line for line in lines if len(line.split())>1 and '# ' not in line]

	exe= ""
	args= ""
	log_foldername= ""
	mapper_infiles= ""
	mapper_outfiles= ""
	reducer_infiles= ""
	reducer_outfiles= ""
	queue_num=-1
	num_queue_line=0

	for line in lines:

		cols = string.split( line )

		if(len(cols)<2):  master_kill_all_slave_jobs_and_exit("len(cols)<2 for line (%s)\nline_number (%d) " %(line,line_number))

		if( cols[0] == "Queue"    ):
			if( len( cols )!= 2 ): master_kill_all_slave_jobs_and_exit("len( cols )!=2 for Queue line (%s)!" %(line))
			num_queue_line+=1
			queue_num = int( cols[1] )
			continue

		####OK, all the other lines should have (cols[1] == "=")####
		if( cols[1] != "=" ): master_kill_all_slave_jobs_and_exit("cols[1] != \"=\" for line (%s)" %(line))

		if( cols[0] ==  "executable" ):
			if( len(cols) != 3 ): master_kill_all_slave_jobs_and_exit("len(cols) != 3 for line (%s)" %(line))
			exe = cols[2]

		if( cols[0] == "#log_foldername" ):
			if( len(cols) != 3 ): master_kill_all_slave_jobs_and_exit("len(cols) != 3 for line (%s)" %(line))
			log_foldername = cols[2]

		if( cols[0] == "arguments"        ): args             = string.join(cols[2:])

		if( cols[0] == "#mapper_outfiles"  ): mapper_outfiles  = string.join(cols[2:])

		if( cols[0] == "#reducer_outfiles" ): reducer_outfiles = string.join(cols[2:])

		if( cols[0] == "#mapper_infiles"   ): mapper_infiles   = string.join(cols[2:]) #NOT CURRENLTY IN USE

		if( cols[0] == "#reducer_infiles"  ): reducer_infiles  = string.join(cols[2:]) #NOT CURRENLTY IN USE

	print "-----------------------------------------------------------------------------------------"
	print "exe=%s" %(exe)
	print "args=%s" %(args)
	print "log_foldername=%s" %(log_foldername)
	print "mapper_infiles=%s" %(mapper_infiles)
	print "mapper_outfiles=%s" %(mapper_outfiles)
	print "reducer_infiles=%s" %(reducer_infiles)
	print "reducer_outfiles=%s" %(reducer_outfiles)
	print "queue_num=%s" %(queue_num)
	print "num_queue_line=%s" %(num_queue_line)
	print "-----------------------------------------------------------------------------------------"

	if(num_queue_line!=1):  master_kill_all_slave_jobs_and_exit("num_queue_line=(%s)!=1" %(num_queue_line))
	if(queue_num==-1): master_kill_all_slave_jobs_and_exit("queue_num==-1")

	reducer_outfile_already_exists = False 
	reducer_outfile_list=reducer_outfiles.split()
	for reducer_outfile in reducer_outfile_list:
		if(exists(reducer_outfile)): 
			reducer_outfile_already_exists = True 
			#master_kill_all_slave_jobs_and_exit("reducer_outfile %s already exist before job submission!" %(reducer_outfile) )


	# Stolen from script for PBS.
	all_job_info_list = []
	#actually_submit_job_info_list = []

	COMPUTER_CLUSTER_NAME=str(os.environ['COMPUTER_CLUSTER_NAME']) #June 15, 2012 TEST.

	if(COMPUTER_CLUSTER_NAME=="LONESTAR-TACC-XSEDE"):
		print "COMPUTER_CLUSTER_NAME==\"LONESTAR-TACC-XSEDE\", extra sleep 0.1 seconds between mkdir of outfile folders"

	if reducer_outfile_already_exists:
		start = queue_num
	else:
		start = 0

	for q in range(start, queue_num + 1 ):

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

			#if(reducer_job): master_kill_all_slave_jobs_and_exit("done_signal for reducer_job %s already exist before job submission!" %(done_signal_filename) )

			for outfile in outfile_list:

				if(outfile[0:7]=="FOLDER:"): #Oct 22, 2011
					outfolder=outfile[7:]
					if(exists(outfolder)==False):	print "WARNING: outfile (%s) doesn't exist" %(outfolder)
				else:
					if(exists(outfile)==False):	  print "WARNING: outfile (%s) doesn't exist" %(outfile)


				print "job %s is already successfully_completed " %(outfile)


		else:
			#######MAKE SURE THAT outfiles and logs files from previously failed runs are deleted####
			already_done=False

			if(exists(out_log_filename)):
				print ("job %s WAS NOT successfully_completed. Removing out_log_filename = %s " %(out_log_filename, out_log_filename))
				master_submit_subprocess( 'rm %s' %(out_log_filename))

			if(exists(err_log_filename)):
				print ("job %s WAS NOT successfully_completed. Removing err_log_filename = %s " %(out_log_filename, err_log_filename))
				master_submit_subprocess( 'rm %s' %(err_log_filename))

			if(exists(job_script_filename)):
				print ("job %s WAS NOT successfully_completed. Removing job_script_filename = %s " %(out_log_filename, job_script_filename))
				master_submit_subprocess( 'rm %s' %(job_script_filename))

			for outfile in outfile_list:

				if(outfile[0:7]=="FOLDER:"): #Oct 22, 2011

					outfolder=outfile[7:]

					if(exists(outfolder)):
						print ("job %s WAS NOT successfully_completed.	 Removing outfolder = %s " %(out_log_filename, outfolder))
						master_submit_subprocess( 'rm -r %s' %(outfolder) )

				else:
					if(exists(outfile)):
						print ("job %s WAS NOT successfully_completed.	 Removing outfile = %s " %(out_log_filename, outfile))
						master_submit_subprocess( 'rm %s' %(outfile) )

			######CREATE THE JOB_SCRIPT FOR THIS JOB #######################################################

			master_safe_mkdir( log_foldername )
			master_safe_mkdir( dirname(out_log_filename) )
			master_safe_mkdir( dirname(err_log_filename) )
			master_safe_mkdir( dirname(job_script_filename) )
			master_safe_mkdir( dirname(done_signal_filename) )

			for outfile in outfile_list:

				if(COMPUTER_CLUSTER_NAME=="LONESTAR-TACC-XSEDE"): sleep(0.1) #June 15, 2012.

				if(outfile[0:7]=="FOLDER:"): #Oct 22, 2011
					outfolder=outfile[7:]
					master_safe_mkdir( dirname(outfolder) )
				else:
					master_safe_mkdir( dirname(outfile) )

			fid = master_safe_open( job_script_filename, 'w' )

			if(reducer_job):
				fid.write('wrapper_script = %s\n' %(get_PYEXE("dagman/DAG_reducer.py")) )
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

	if( (DAG_submit_file_.count('cluster')>0) or (DAG_submit_file_.count('CLUSTER')>0) ):
		Is_clusterer_job=True
		print "%s is a clusterer job!" %(DAG_submit_file_)
	else:
		Is_clusterer_job=False
		print "%s is NOT a clusterer job!" %(DAG_submit_file_)

	all_job_info_list=submit_jobs_to_slave( all_job_info_list, Is_clusterer_job )



#	if(queue_num!=0):
#	else:
#		if(len(all_job_info_list)!=0): master_kill_all_slave_jobs_and_exit("queue_num==0 but len(all_job_info_list=0")
#		#consistency check
#		print "Queue_num=0 for JOB %s, no actual job submission!" %(DAG_submit_file_)


	print "A total of %f seconds is used to submit the JOB: %s " %( time.time()- DAG_submit_start_time, DAG_submit_file_)

	return all_job_info_list


