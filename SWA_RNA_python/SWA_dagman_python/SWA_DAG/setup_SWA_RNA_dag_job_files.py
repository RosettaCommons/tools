#!/usr/bin/env python

######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options

start_argv=copy.deepcopy(argv)
######################################################################


#setup_SWA_RNA_dag_job_files.py -s template.pdb -fasta fasta -sample_res 3-8 -single_stranded_loop_mode True -local_demo True -native_pdb native.pdb

#setup_SWA_RNA_dag_job_files.py -s small_template.pdb -fasta small_fasta  -sample_res 3-8 -single_stranded_loop_mode True -local_demo True -native_pdb  small_native.pdb 

####################################################################################

local_demo= parse_options( argv, "local_demo", "False", Verbose=False )

command= "create_dag_job_files.py %s > LOG_create_dag_job_files.out" %(list_to_string(argv[1:]))

submit_subprocess(command)

submit_subprocess("python README_SETUP.py")

if(local_demo):
	local_demo_sampler_command	 	 =get_PYEXE("misc/create_local_dag.py")
	local_demo_sampler_command	 	+=" -dag_file CONDOR/SAMPLER/REGION_0_1_START_FROM_REGION_0_0.condor "
	local_demo_sampler_command	 	+=" -folder_layer 0 -local_dag_file SAMPLER_DEMO_COMMAND -minimizer_rename_tag True -VERBOSE False -output_pdb False > LOG_create_local_demo_sampler_command.txt "
	submit_subprocess(local_demo_sampler_command)

	local_demo_clusterer_command	 =get_PYEXE("misc/create_local_dag.py")
	local_demo_clusterer_command	+=" -dag_file  CONDOR/CLUSTERER/REGION_0_1_cluster.condor -clusterer_silent_file_in region_0_1_sample.out "
	local_demo_clusterer_command	+=" -folder_layer 0 -local_dag_file CLUSTERER_DEMO_COMMAND -minimizer_rename_tag True -VERBOSE False -output_pdb False > LOG_create_local_demo_clusterer_command.txt "
	submit_subprocess(local_demo_clusterer_command)

	README_SETUP = open( "LOCAL_DEMO", 'w')
	README_SETUP.write( '\nsource SAMPLER_DEMO_COMMAND\n' )
	README_SETUP.write( '\nsource CLUSTERER_DEMO_COMMAND\n' )
	README_SETUP.close()

####################################################################################
print "----------------------------------------------------------------------------------------------------------------------------"
print "Successfully RAN %s" %( list_to_string(start_argv) )
print "----------------------------------------------------------------------------------------------------------------------------"


