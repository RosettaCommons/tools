#!/usr/bin/env python

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
import os
from os.path import basename, dirname, exists, expanduser
from time import sleep
from sets import Set
import copy
######################################################################

from FARFAR_dag_util import *

######################################################################
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.utility.DAGMAN_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################

HOMEDIR = expanduser('~')

#########################################################################################################


PRE_PROCESS_SCRIPT = get_PYEXE("FARFAR_DAG/FARFAR_pre_process.py")
FARFAR_FILTERER_EXE= get_PYEXE("FARFAR_DAG/FARFAR_filterer.py")
FARFAR_FINAL_FILTER_EXE=get_PYEXE("FARFAR_DAG/FARFAR_final_filter_dag.py")
SWA_CLUSTER_EXE=get_PYEXE("misc/SWA_cluster.py")
GENERIC_EMPTY_REDUCER_SCRIPT= get_PYEXE("dagman/generic_empty_reducer.py")



#DAG_FILE_FOLDER=abspath("CONDOR/")
DAG_FILE_FOLDER="CONDOR/" #Aug 24, 2011..new get_DAG_log_foldername() is sensitive to dirname!

COMMON_ARGS_FILE=get_FARFAR_common_args_file()

CLUSTER_ARGS_FILE=get_FARFAR_cluster_args_file()

ROSETTA_COMMAND_FILE=get_FARFAR_ROSETTA_COMMAND_file()

if(exists(COMMON_ARGS_FILE)==False): error_exit_with_message("COMMON_ARGS_FILE (%s) doesn't exist!" %(COMMON_ARGS_FILE) )

if(exists(CLUSTER_ARGS_FILE)==False): error_exit_with_message("CLUSTER_ARGS_FILE (%s) doesn't exist!" %(CLUSTER_ARGS_FILE) )

if(exists(ROSETTA_COMMAND_FILE)==False): error_exit_with_message("ROSETTA_COMMAND_FILE (%s) doesn't exist!" %(ROSETTA_COMMAND_FILE))

#########################################################################################################

def setup_FINAL_CLUSTER_job(fid_dag, FINAL_CLUSTER_job_tag, FINAL_FILTER_tag, final_filtered_file):

	extra_cluster_args_string=safe_readlines( CLUSTER_ARGS_FILE )[0][:-1]

	FINAL_CLUSTER_filename= '%s/%s.condor' %(DAG_FILE_FOLDER, FINAL_CLUSTER_job_tag)
	clustered_outfile="clustered_%s" %(final_filtered_file)

	if(exists(clustered_outfile)==False):
	
		cluster_arguments=  " -silent_file %s " %(final_filtered_file)
		cluster_arguments+= " -output_filename %s " %(clustered_outfile)
		cluster_arguments+= " -cluster_rmsd 0.7 "
		cluster_arguments+= " -suite_cluster_rmsd 1.0 "
		cluster_arguments+= " -suite_cluster_OPTION True "
		cluster_arguments+= " -distinguish_pucker true "
		cluster_arguments+= " -extract_pdb False "
		cluster_arguments+= " -clusterer_rename_tags true "
		cluster_arguments+= " -add_lead_zero_to_tag false "
		cluster_arguments+= " -num_pose_kept 10000 "
		cluster_arguments+= " -VERBOSE false "
		cluster_arguments+= " %s " %(extra_cluster_args_string)
		cluster_arguments+= " -common_args %s " %(COMMON_ARGS_FILE)
		cluster_arguments+= " -redirect_out_log False "

		make_dag_job_submit_file( FINAL_CLUSTER_filename, SWA_CLUSTER_EXE, cluster_arguments, '', clustered_outfile, '', '')

		fid_dag.write('\nJOB %s %s\n' % ( FINAL_CLUSTER_job_tag, FINAL_CLUSTER_filename) )
		if(exists(final_filtered_file)==False): fid_dag.write('PARENT %s CHILD %s\n' %(FINAL_FILTER_tag, FINAL_CLUSTER_job_tag) )


	return clustered_outfile

