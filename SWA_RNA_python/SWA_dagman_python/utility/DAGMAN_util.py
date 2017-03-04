#!/usr/bin/env python

###########################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options

###########################

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


#############################################################################################################
def create_generic_README_SUB(num_slave_nodes):

	if(isinstance( num_slave_nodes, int )==False):
		print "PROBLEM num_slave_nodes=", num_slave_nodes
		error_exit_with_message("num_slave_nodes object is not an int!")

	master_script_EXE=get_PYEXE("dagman/DAG_continuous.py")

	README_SUB = open( "README_SUB.py", 'w')

	README_SUB.write( '#!/usr/bin/env python\n' )
	README_SUB.write( 'from os import system\n' )
	README_SUB.write( 'from os.path import exists\n' )
	README_SUB.write( 'from os.path import abspath\n' )

	README_SUB.write( 'import string\n\n' )
	README_SUB.write( 'if( exists("master_log.out") ): system("rm master_log.out")\n' )
	README_SUB.write( 'if( exists("master_log.err") ): system("rm master_log.err")\n' )
	README_SUB.write( 'if( exists("SLAVE_JOBS/") ): system("rm -rf SLAVE_JOBS/")\n\n' )
	README_SUB.write( 'master_tag = abspath( "MASTER" ).replace("/","_")\n\n')

	command_1=	"bsub -W 144:0 -o master_log.out -e master_log.err "
	command_2= " %s  -j %d rna_build.dag" %(master_script_EXE, num_slave_nodes)

	README_SUB.write( 'command_act=' + "\"" + command_1 + "\" + " + '(\" -J %s \" %(master_tag))' " + \"" + command_2 + "\"" + '\n')

	#README_SUB.write('print "command_act=", command_act')
	README_SUB.write( 'system(command_act)\n' )
	README_SUB.close()

#############################################################################################################
def get_chemical_shift_args_file():

	return "CHEMICAL_SHIFT_ARGS.txt"

#############################################################################################################
def parse_chemical_shift_args( chemical_shift_args_file ):

	chemical_shift_arguments=" "

	chemical_shift_args=safe_readlines( chemical_shift_args_file )[0].split()

	H5_prime_mode= parse_options( chemical_shift_args, "H5_prime_mode", "")

	if(H5_prime_mode==""): error_exit_with_message("H5_prime_mode==\"\"!!")

	chemical_shift_arguments+=" -H5_prime_mode %s " %(H5_prime_mode)

	return chemical_shift_arguments

#############################################################################################################
def setup_chemical_shift_args( argv ):

	chemical_shift_H5_prime_mode = parse_options( argv, "chemical_shift_H5_prime_mode", "LEAST_SQUARE_IGNORE_DUPLICATE" ) #Before Nov 19, 2011. Default is LEAST_SQUARE

	chemical_shift_args_file=get_chemical_shift_args_file()

	if(exists(chemical_shift_args_file)):
		print "Warning chemical_shift_args_file (%s) already exist......removing!" %(chemical_shift_args_file)
		submit_subprocess("rm %s" %(chemical_shift_args_file) )
		#error_exit_with_message("chemical_shift_args_file (%s) already exist!" %(chemical_shift_args_file))

	fid = open( chemical_shift_args_file, 'w' )
	fid.write(' -H5_prime_mode %s ' %(chemical_shift_H5_prime_mode) )
	fid.close()

#############################################################################################################

def get_condor_submit_file( job_tag , TYPE=""):

	if(TYPE==""):
		condor_submit_file='CONDOR/%s.condor' %(job_tag)
	else:
		condor_submit_file='CONDOR/%s/%s.condor' %(TYPE,job_tag)


	return condor_submit_file

##############################################################################################################


def get_DAG_log_foldername(dag_job_filename):

	#log_foldername = 'CONDOR/' + basename(dag_job_filename).replace('.condor','') #A example is the folder REGION_0_1_START_FROM_REGION_0_0/
	log_foldername = dag_job_filename[:-7] #Change to this on July 20, 2011

	return log_foldername


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

	# New: needed for 'regular' condor DAGMAN submission
	job_info = get_job_info( log_foldername, '$(Process)' )
	fid.write('output = %s\n' % job_info[ 'out_log_filename'] )
	fid.write('error = %s\n' % job_info[ 'err_log_filename' ] )
	fid.write('log = %s\n' % 'condor.log' )
	fid.write('Notification = error\n' )

	fid.write('Queue %d\n' % 1 )

	# Note from rhiju -- adding hashtag at beginning permits these files to be used by condor DAGman (condor_submit_dag)
	fid.write('\n')
	fid.write('################################################################\n')
	fid.write('# Below lines are commented so that they don''t disrupt runs \n' )
	fid.write('# that are queued by condor DAGman. But they are in-use by other scripts. \n' )
	fid.write('################################################################\n')
	fid.write('#mapper_infiles = %s\n' % mapper_infiles) #currently used only in SWA2
	fid.write('#mapper_outfiles = %s\n' % mapper_outfiles)
	fid.write('#reducer_infiles = %s\n' % reducer_infiles) #currently used only in SWA2
	fid.write('#reducer_outfiles = %s\n' % reducer_outfiles)
	fid.write('#log_foldername = %s\n' % log_foldername)
	fid.close()


##########This is used in the pre_process scripts#################################
def update_CONDOR_file_with_actual_job_queue_num(condor_submit_file, N_JOBS): #Updated on FEB 08, 2012

	print "%d JOBS for condor_submit_file %s" %(N_JOBS, condor_submit_file)

	if(condor_submit_file==""): error_exit_with_message("condor_submit_file==\"\"!")

	if(exists(condor_submit_file)==False): error_exit_with_message("condor_submit_file (%s) doesn't exist!" %(condor_submit_file) )

	# Go through condor submission file and update number of jobs...
	lines = open( condor_submit_file ).readlines()

	submit_subprocess("rm %s" %(condor_submit_file) )

	fid = open( condor_submit_file, 'w' )

	num_queue_line=0

	for line in lines:

		cols=line.split()

		if( len(cols)<2 ): continue # error_exit_with_message("len(cols)<2 for line (%s)" %(line))

		if(cols[0] == 'Queue'):

			num_queue_line+=1

			if( len(cols)!=2 ): error_exit_with_message("len(cols)!=2 for Queue line=(%s)" %(line))

			fid.write( 'Queue %d\n' % N_JOBS )

		else:

			fid.write( line )

	fid.close()

	if(num_queue_line!=1): error_exit_with_message("num_queue_line=(%s)!=1" %(num_queue_line))

##############################################################################################################################
def make_common_args_file( specific_common_args, input_filename ):

	foldername='COMMON_ARGS/'
	filename= foldername + input_filename

	if (exists( foldername )==False):
		system( 'mkdir -p ' + foldername )

	if (exists( filename )):
		system( 'rm ' + filename ) #Actually not necessary

	fid = open( filename, 'w' )
	fid.write( specific_common_args )
	fid.close()


################################################################################################################################################
def setup_final_rebuild_bulge_dag_job_file(fid_dag, input_silent_file, rebuild_bulge_job_tag, native_pdb, \
																				 common_args_file, main_algorithm, BMRB_chemical_shift_file):

	##########################################################################################################################################
	DAG_REBUILD_BULGE_PRE_PROCESS=get_PYPATH("SWA_DAG/DAG_rebuild_bulge_pre_process.py")
	DAG_REBUILD_BULGE_MAPPER= get_PYPATH("SWA_DAG/DAG_rebuild_bulge.py")
	DAG_REBUILD_BULGE_REDUCER= get_PYPATH("SWA_DAG/DAG_rebuild_bulge_reducer.py")

	if( PATH_exists( DAG_REBUILD_BULGE_PRE_PROCESS ) == False ):  error_exit_with_message("DAG_REBUILD_BULGE_PRE_PROCESS (%s) doesn't exist!" %(DAG_REBUILD_BULGE_PRE_PROCESS) )
	if( PATH_exists( DAG_REBUILD_BULGE_MAPPER ) == False ):  error_exit_with_message("DAG_REBUILD_BULGE_MAPPER (%s) doesn't exist!" %(DAG_REBUILD_BULGE_MAPPER) )
	if( PATH_exists( DAG_REBUILD_BULGE_REDUCER ) == False ):  error_exit_with_message("DAG_REBUILD_BULGE_REDUCER (%s) doesn't exist!" %(DAG_REBUILD_BULGE_REDUCER) )
	##########################################################################################################################################

	if(exists(common_args_file)==False): error_exit_with_message("common_args_file (%s) doesn't exist!" %(common_args_file) )


	condor_submit_file=get_condor_submit_file(rebuild_bulge_job_tag, "REBUILD_BULGE")

	log_foldername = get_DAG_log_foldername(condor_submit_file)

	reducer_job_info = get_job_info(log_foldername, 'reducer')

	############################################################################################################################################
	if(rebuild_bulge_job_tag.count("WITH_CHEM_SHIFT")==0):

		outfolder="%s/STANDARD/" %(rebuild_bulge_job_tag)

		if(BMRB_chemical_shift_file!=""): error_exit_with_message("rebuild_bulge_job_tag==\"FINAL_REBUILD_BULGE\" but BMRB_chemical_shift_file!=\"\"")

	elif(rebuild_bulge_job_tag.count("WITH_CHEM_SHIFT")==1):

		if(rebuild_bulge_job_tag[-16:]!="_WITH_CHEM_SHIFT"):
			error_exit_with_message("rebuild_bulge_job_tag[-16:]!=\"_WITH_CHEM_SHIFT\" for rebuild_bulge_job_tag (%s)" %(rebuild_bulge_job_tag))

		outfolder="%s/WITH_CHEM_SHIFT/" %(rebuild_bulge_job_tag[:-16])

		if(BMRB_chemical_shift_file==""): error_exit_with_message("rebuild_bulge_job_tag==\"FINAL_REBUILD_BULGE_WITH_CHEM_SHIFT\" but BMRB_chemical_shift_file==\"\"")

		if(exists(BMRB_chemical_shift_file)==False): error_exit_with_message("BMRMB_chemical_shift_file (%s) doesn't exist!" %(BMRB_chemical_shift_file) )

	else:
		error_exit_with_message("Invalid rebuild_bulge_job_tag (%s)" %(rebuild_bulge_job_tag))
	############################################################################################################################################

	submit_subprocess("mkdir -p %s" %(outfolder))

	reducer_outfile="%s/rebuild_bulge_%s" %(outfolder, basename(input_silent_file))

	######OK DAG is done! Important, cannot submit DAG that already done, since DAG_continuous.py will raise error if reducer_done_signal_filename already exist!####
	if( exists(reducer_job_info['done_signal_filename']) ):
		if(exists(reducer_outfile)==False): error_exit_with_message("reducer_outfile (%s) doesn't exist!" %(reducer_outfile) )
		return False

	#####OK if reducer_outfile exist then assume that job is done#####Note that this behavior differ from job-types to job-types####
	#####See create_sampled_virt_ribose_silent_file() for the oppositve behavior####################################################
	if( exists(reducer_outfile) ):
		print "reducer_outfile (%s) already exist! Will not submit job (%s) " %(reducer_outfile, rebuild_bulge_job_tag)
		return False

	fid_dag.write('\nJOB %s %s\n' % ( rebuild_bulge_job_tag, condor_submit_file) )

	##################PRE_PROCESS#################

	pre_process_command="%s -silent_file %s -condor_submit_file %s " %(DAG_REBUILD_BULGE_PRE_PROCESS, input_silent_file, condor_submit_file)

	fid_dag.write('SCRIPT PRE %s %s\n' % (rebuild_bulge_job_tag, pre_process_command) )

	##################REDUCER#################

	reducer_command="%s -outfolder %s -reducer_outfile %s -condor_submit_file %s " %(DAG_REBUILD_BULGE_REDUCER, outfolder, basename(reducer_outfile), condor_submit_file)

	fid_dag.write('SCRIPT POST %s %s \n' % (rebuild_bulge_job_tag, reducer_command ) )

	queue_tag = 'S_$(Process)'  ###umm this specific form is assumed in SWA_sampling_pre_process.py

	##################MAPPER#######################

	mapper_outfolder="%s/%s/" %(outfolder, queue_tag)

	mapper_outfile="%s/mapper_rebuild_bulge.out" %(mapper_outfolder)

	mapper_arguments=" -start_silent_file %s -common_args_file %s " %(input_silent_file, common_args_file)
	mapper_arguments+=" -queue_ID $(Process) -outfolder %s -out:file:silent %s " %(mapper_outfolder, mapper_outfile)
	mapper_arguments+=" -native_pdb %s -main_algorithm %s " %(native_pdb, main_algorithm)

	if(BMRB_chemical_shift_file!=""):
		mapper_arguments+=" -BMRB_chemical_shift_file %s " %(BMRB_chemical_shift_file)

		mapper_arguments+=" -chemical_shift_args_file %s " %( get_chemical_shift_args_file() )

	#########Should always include the mapper_outfile before the mapper_outfolder???#####################
	make_dag_job_submit_file( condor_submit_file, DAG_REBUILD_BULGE_MAPPER, mapper_arguments, '', "%s FOLDER:%s " %(mapper_outfile, mapper_outfolder) , '', reducer_outfile)

	return True
################################################################################################################################################

