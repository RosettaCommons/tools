#!/usr/bin/python

######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options

start_argv=copy.deepcopy(argv)
######################################################################

###This is the RELEASE SCRIPT###

#setup_SWA_RNA_dag_job_files.py -s template.pdb -fasta C7_2_target.fasta -sample_res 11 12 13 -nstruct 1000 -num_slave_nodes 250 -single_stranded_loop_mode True

#setup_SWA_RNA_dag_job_files.py -s template.pdb -sample_res 3-8 -nstruct 1000 -num_slave_nodes 500 -single_stranded_loop_mode True -native_pdb native.pdb -native_virtual_res 5

#-nstruct 1000 -num_slave_nodes 500 -native_virtual_res 5 

#setup_SWA_RNA_dag_job_files.py -s template.pdb -sample_res 3-8 -single_stranded_loop_mode True -native_pdb native.pdb -force_field_file rna/rna_loop_hires_04092010.wts -rna_torsion_potential_folder ps_03242010/ -sample_virt_ribose_in_sep_DAG True -apply_VDW_rep_delete_matching_res True -clusterer_optimize_memory_usage true  -local_demo True  -fasta fasta 


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
	local_demo_sampler_command	 	+=" -folder_layer 0 -local_dag_file SAMPLER_DEMO_COMMAND -minimizer_rename_tags True -VERBOSE False -output_pdb False > LOG_create_local_demo_sampler_command.txt "
	submit_subprocess(local_demo_sampler_command)

	local_demo_clusterer_command	 =get_PYEXE("misc/create_local_dag.py")
	local_demo_clusterer_command	+=" -dag_file  CONDOR/CLUSTERER/REGION_0_1_cluster.condor -clusterer_silent_file_in region_0_1_sample.out "
	local_demo_clusterer_command	+=" -folder_layer 0 -local_dag_file CLUSTERER_DEMO_COMMAND -minimizer_rename_tags True -VERBOSE False -output_pdb False > LOG_create_local_demo_clusterer_command.txt "
	submit_subprocess(local_demo_clusterer_command)


####################################################################################
print "----------------------------------------------------------------------------------------------------------------------------"
print "Successfully RAN %s" %( list_to_string(start_argv) )
print "----------------------------------------------------------------------------------------------------------------------------"



































######NOTE TO SELF#####
##DO TO:
###0. Quick test mode (run the first SWA building step i.e. region 0_1 | meant for testing on a local computer)
###1. Handle the Possibility of no native_pdb! (except fasta instead!)
###2. RENAME native_virtual_res to native_bulge_res (need to make corresponding changes in create_benchmark_job_files and setup_FARFAR_setup_files.py
###3. call the Rosetta_ready and renumber_pdb /replace chain in place scripts.


#These options are defualt in (is_release_mode() and single_stranded_loop_mode)
#" -force_field_file rna/rna_loop_hires_04092010.wts -rna_torsion_potential_folder ps_03242010/ -sample_virt_ribose_in_sep_DAG True -apply_VDW_rep_delete_matching_res True "

##RENAME:
#expand_radius_100_1s72_RNA_A_35_40_3_8.pdb -->native.pdb
#no_loop_expand_radius_100_1s72_RNA_A_35_40_3_8.pdb -->template.pdb
#no_loop_expand_radius_300_1s72_RNA_A_35_40.pdb 	-->peripheral.pdb

#Extra (optional) flags:
#-native_pdb (ex. native.pdb)
#-native_virtual_res (seq_num) #RMSD calculation purposes!
#RMSD_MODE?

#-VDW_rep_screen_info (ex. -VDW_rep_screen_info peripheral.pdb 	12-187 1-47)
#-excise_segment_torsion_DB (ex. 34-41 (Extra FARFAR option!)) #Correspond to seq_num 34-41 of the 1jj2.torsions!
#-fragment_match_type MATCH_YR (Extra FARFAR option!..don't include this!)

