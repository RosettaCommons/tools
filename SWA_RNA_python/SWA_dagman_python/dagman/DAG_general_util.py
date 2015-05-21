#!/usr/bin/env python

###########################

from SWA_dagman_python.utility.SWA_util import *

###########################

#######General utility function related to generic DAG job submission#############

#######For DAG stuff specific to SWA_RNA_DAG or FARFAR_RNA_DAG, should put in utility/RNA_dagman_util.py instead######


#############################################################################################################
def create_generic_README_SUB(num_slave_nodes, wall_time, readme_sub_filename):

	if(isinstance( num_slave_nodes, int )==False): 
		print "PROBLEM num_slave_nodes=", num_slave_nodes
		error_exit_with_message("num_slave_nodes object is not an int!")

	if(isinstance( wall_time, int )==False): 
		print "PROBLEM wall_time=", wall_time
		error_exit_with_message("wall_time object is not an int!")

	submit_DAG_job_EXE=get_PYEXE("dagman/submit_DAG_job.py")

	README_SUB = open( readme_sub_filename, 'w')

	README_SUB.write( '#!/usr/bin/env python\n' )
	README_SUB.write( 'from os import system\n' )
	#README_SUB.write( 'from os.path import exists\n' )
	#README_SUB.write( 'from os.path import abspath\n' )

	README_SUB.write( 'import string\n\n' )

	#job_name, job_dir_name, outfile, errfile, job_script, walltime=168, 	
	#Before April 27, 2012 used to be -walltime 144

	README_SUB.write( 'command=  "%s " \n\n' %(submit_DAG_job_EXE))

	README_SUB.write( 'command+= "-dagman_file rna_build.dag " \n\n')

	README_SUB.write( 'command+= "-master_memory_reserve 2048 " \n\n')

	README_SUB.write( 'command+= "-master_wall_time %d " \n\n' %(wall_time))

	README_SUB.write( 'command+= "-num_slave_nodes %d " \n\n' %(num_slave_nodes))

	README_SUB.write( 'command+= "> LOG_submit_DAG_job.txt " \n\n')

	README_SUB.write( 'system(command)\n' )

	README_SUB.close()


######################################################################
def	get_job_info(log_foldername, job_name):

	job_info = {}

	prefix=""
	if(job_name=='reducer' or job_name=='0'):
		prefix='KEEP_LOG_FILE/'

	job_info['out_log_filename'] = prefix + log_foldername + '/outfile/' + job_name + '.out'
	job_info['err_log_filename'] = prefix + log_foldername + '/errfile/' + job_name + '.err'
	job_info['job_script_filename'] = prefix + log_foldername + '/job_script/' + job_name + '.qsub'
	job_info['done_signal_filename'] = 'DONE/' + log_foldername + '/' + job_name + '.txt' 

	return job_info

######################################################################
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


##############################################################################################################
def get_DAG_log_foldername(dag_job_filename):

	#log_foldername = 'CONDOR/' + basename(dag_job_filename).replace('.condor','') #A example is the folder REGION_0_1_START_FROM_REGION_0_0/
	log_foldername = dag_job_filename[:-7] #Change to this on July 20, 2011

	return log_foldername

#############################################################################################################
def get_dag_job_submit_file( job_tag , TYPE=""): 

	if(TYPE==""):
		condor_submit_file='CONDOR/%s.condor' %(job_tag)
	else:
		condor_submit_file='CONDOR/%s/%s.condor' %(TYPE,job_tag)


	return condor_submit_file

##############################################################################################################

def make_dag_job_submit_file( dag_job_filename, EXE, arguments, mapper_infiles, mapper_outfiles, reducer_infiles, reducer_outfiles):

	if(dag_job_filename[-7:]!='.condor'): error_exit_with_message("dag_job_filename (%s) need to have a '.condor' extension!" %(dag_job_filename) ) 

	if exists( dag_job_filename ):
		print 'rm %s' %(dag_job_filename) 
		system( 'rm %s' %(dag_job_filename) ) #remove a old version of the condor_submit_file if it happens to exist		

	log_foldername = get_DAG_log_foldername(dag_job_filename)	

	if not exists( log_foldername ): system( 'mkdir -p ' + log_foldername )


	fid = open( dag_job_filename, 'w' )

	fid.write('executable = %s\n' % EXE )
	fid.write('arguments = %s\n' % arguments)
	fid.write('mapper_infiles = %s\n' % mapper_infiles) #currently used only in SWA2
	fid.write('mapper_outfiles = %s\n' % mapper_outfiles) 
	fid.write('reducer_infiles = %s\n' % reducer_infiles) #currently used only in SWA2
	fid.write('reducer_outfiles = %s\n' % reducer_outfiles) 
	fid.write('log_foldername = %s\n' % log_foldername)
	fid.write('Queue %d\n' % 1 )
	fid.close()
