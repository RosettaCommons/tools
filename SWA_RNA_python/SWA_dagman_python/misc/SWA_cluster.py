#!/usr/bin/env python

from os.path import exists,dirname,basename,expanduser,abspath
from sys import exit, argv
import string
from time import sleep
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import *
######################################################################


print_title_text("Enter: " + list_to_string(copy.deepcopy(argv)))
print
print

#

#SWA_cluster.py -silent_file    FINAL_filtered_energy.out   -cluster_rmsd 0.0 -skip_clustering True -extract_pdb False  -clusterer_rename_tags true -add_lead_zero_to_tag false 

#SWA_cluster.py -silent_file     replace_scoreline_WITH_SHIFT_STATS_rebuild_Jan_13_1JJ2.out     -common_args FARFAR_common_args.txt    -cluster_rmsd 0.0 -skip_clustering True -extract_pdb False -native_pdb 1JJ2_S_6983.pdb      -recreate_silent_struct True -clusterer_rename_tags false

#SWA_cluster.py -silent_file  silent_file.out -add_lead_zero_to_tag false -cluster_rmsd 0.0

#SWA_cluster.py -silent_file  WITH_SHIFT_STATS_Jan_25_RNA09_Jan_28_1JJ2.out -common_args FARFAR_common_args.txt  -cluster_rmsd 2.0 -suite_cluster_rmsd 2.5 -suite_cluster_OPTION True -full_length_loop_rmsd_clustering True -num_pose_kept 50

#SWA_cluster.py -silent_file   recal_rmsd_NO_SHIFT_STATS_combine.out  -common_args MANUAL_common_args.txt -cluster_rmsd 2.0 -suite_cluster_rmsd 2.5 -suite_cluster_OPTION True -full_length_loop_rmsd_clustering True -num_pose_kept 50

#SWA_cluster.py -silent_file WITH_SHIFT_STATS_region_FINAL.out -common_args common_args_region_FINAL.out   -cluster_rmsd 2.0 -suite_cluster_rmsd 2.5 -suite_cluster_OPTION True -full_length_loop_rmsd_clustering True  -num_pose_kept 50

#SWA_cluster.py -silent_file WITH_IMINO_PENALTY.out WITHOUT_IMINO_PENALTY.out  -cluster_rmsd 0.0  -skip_clustering True -extract_pdb False -clusterer_rename_tags false

#SWA_cluster.py -silent_file recal_full_WITH_SHIFT_STATS_region_FINAL.out -use_best_neighboring_shift_RMSD True -common_args common_args_region_FINAL.out  -cluster_rmsd 1.0 -suite_cluster_rmsd 1.5 -suite_cluster_OPTION True -extract_pdb False -skip_clustering True

#SWA_cluster.py -silent_file shift_RMSD_max_n_struct_1000_WITH_SHIFT_STATS_combine.out -use_best_neighboring_shift_RMSD True -common_args ../../MANUAL_common_args.txt  -cluster_rmsd 0.7 -suite_cluster_rmsd 1.0 -suite_cluster_OPTION True -extract_pdb False -skip_clustering True

#SWA_cluster.py -silent_file NEW_O_loop_rmsd_abs_score_cut_150_WITH_SHIFT_STATS_combine.out -use_best_neighboring_shift_RMSD True -common_args ../../MANUAL_common_args.txt  -cluster_rmsd 1.5  -suite_cluster_OPTION True -extract_pdb False -skip_clustering True

#SWA_cluster.py -silent_file REWEIGHT_shift_score_3.00_WITH_SHIFT_STATS_combine.out -common_args MANUAL_common_args.txt  -cluster_rmsd 2.0 -suite_cluster_rmsd 2.5 -suite_cluster_OPTION True -full_length_loop_rmsd_clustering True 

#SWA_cluster.py -silent_file shift_RMSD_abs_score_cut_30_WITH_SHIFT_STATS_combine.out -common_args ../common_args_region_FINAL.out -cluster_rmsd 2.0 -suite_cluster_rmsd 2.5 -suite_cluster_OPTION True -full_length_loop_rmsd_clustering True 

#SWA_cluster.py -silent_file NEW_O_loop_rmsd_abs_score_cut_150_shift_RMSD_abs_score_cut_25_WITH_SHIFT_STATS_combine.out -common_args MANUAL_common_args.txt -cluster_rmsd 2.0 -suite_cluster_rmsd 2.5 -suite_cluster_OPTION True 

#SWA_cluster.py -silent_file recal_rmsd_WITH_SHIFT_STATS_combine.out -common_args MANUAL_common_args.txt  -cluster_rmsd 0.7 -suite_cluster_rmsd 1.0 -suite_cluster_OPTION True -extract_pdb False -clusterer_rename_tags false -num_pose_kept 500000

#SWA_cluster.py -silent_file WITH_SHIFT_STATS_combine.out -common_args ../common_args_region_FINAL.out  -cluster_rmsd 0.0 -skip_clustering True -extract_pdb False -native_pdb SWA_S_759.pdb  -recreate_silent_struct True -clusterer_rename_tags false

#SWA_cluster.py -silent_file WITH_SHIFT_STATS_combine.out -common_args common_args_region_FINAL.out   -cluster_rmsd 0.0 -skip_clustering True -extract_pdb False -native_pdb M_1.pdb   -recreate_silent_struct True -clusterer_rename_tags false

#SWA_cluster.py -silent_file WITH_SHIFT_STATS_FARFAR.out -common_args ../../FARFAR_common_args.txt  -cluster_rmsd 0.0 -skip_clustering True -extract_pdb False -native_pdb M_1.pdb -recreate_silent_struct True -clusterer_rename_tags false


#SWA_cluster.py -silent_file WITH_SHIFT_STATS_rebuild_bulge_clustered_FINAL_filtered_energy.out -common_args MANUAL_common_args.txt  -cluster_rmsd 0.0 -skip_clustering True -extract_pdb False -native_pdb 3X2_loop_2aht.pdb -recreate_silent_struct True -clusterer_rename_tags false

#SWA_cluster.py -silent_file shift_RMSD_abs_score_cut_25_WITH_SHIFT_STATS_combine.out -common_args MANUAL_common_args.txt -cluster_rmsd 2.0 -suite_cluster_rmsd 2.5 -suite_cluster_OPTION True 

#SWA_cluster.py -silent_file WITH_SHIFT_STATS_rebuild_bulge_region_FINAL.out -recreate_silent_struct True -skip_clustering True -native_pdb HUMAN_S_0.pdb -common_args ../../COMMON_ARGS/common_args_region_FINAL.out  -extract_pdb False  -clusterer_rename_tags false

#SWA_cluster.py -silent_file WITH_SHIFT_STATS_combine.out  -cluster_rmsd 2.0 -suite_cluster_rmsd 2.5 -suite_cluster_OPTION True -common_args common_args_region_FINAL.out 


#SWA_cluster.py -silent_file start_from_region_0_3_sample_filtered.out -cluster_rmsd 0.7 -suite_cluster_rmsd 1.0 -suite_cluster_OPTION True -extract_pdb False -num_pose_kept 5000 -common_args ../COMMON_ARGS/common_args_region_FINAL.out

#SWA_cluster.py -silent_file region_FINAL.out  -cluster_rmsd 0.0 -skip_clustering True -extract_pdb False   -common_args ../COMMON_ARGS/common_args_region_FINAL.out -native_pdb short_loop_1rng.pdb -recreate_silent_struct True 


#SWA_cluster.py -silent_file recal_rmsd_FORCE_SYN_CHI_FARFAR_filter_energy.out recal_rmsd_NO_CONSTRAINT_FARFAR_filter_energy.out  -cluster_rmsd 0.0 -skip_clustering True -extract_pdb False   

#
#
#

#SWA_cluster.py -silent_file REWEIGHT_fa_stack_4.00_rna_bulge_2.22_Oct_7_region_FINAL.out region_FINAL.out -cluster_rmsd 0.0 -skip_clustering True -extract_pdb False -add_lead_zero_to_tag false



#SWA_cluster.py -silent_file cluster_200_shift_RMSD_max_n_struct_377_WITH_SHIFT_RMSD_combine_region_FINAL.out -common_args FILTER_common_args.out -perform_filters True -skip_clustering True -extract_pdb False -clusterer_rename_tags false


#SWA_cluster.py -silent_file FINAL_filtered_energy.out -cluster_rmsd 0.7 -suite_cluster_rmsd 1.0 -suite_cluster_OPTION True -extract_pdb False -num_pose_kept 2000 -common_args common_args_region_FILTER.out

#SWA_cluster.py -silent_file FINAL_filtered_energy.out -cluster_rmsd 1.5 -suite_cluster_rmsd 2.5 -suite_cluster_OPTION True -extract_pdb False -num_pose_kept 2000 -common_args common_args_region_FILTER.out


#SWA_cluster.py -silent_file FINAL_filtered_energy.out -cluster_rmsd 2.0 -extract_pdb False -num_pose_kept 2000 -common_args common_args_region_FILTER.out


#SWA_cluster.py -silent_file FINAL_filtered_energy.out -cluster_rmsd 0.75 -num_pose_kept 2000 


#-suite_cluster_radius 1.0  -loop_cluster_radius 0.7


#SWA_cluster.py -silent_file DAG_ID_0_filtered_energy.out -cluster_rmsd 0.0 -skip_clustering True -extract_pdb False

#SWA_cluster.py -silent_file region_FINAL.out   -common_args COMMON_ARGS/common_args_region_FINAL.out -native_pdb June_05_S_0.pdb -recreate_silent_struct True -skip_clustering True -extract_pdb False


#SWA_cluster.py -silent_file ROSETTA_filter_region_FINAL.out  -common_args common_args_region_FILTER.out -cluster_rmsd 4.0 -suite_cluster_rmsd 4.0 -suite_cluster_OPTION True 

#SWA_cluster.py -silent_file FINAL_filtered_energy.out -common_args common_args_region_FILTER.out -perform_filters True -skip_clustering True -extract_pdb False


#SWA_cluster.py -silent_file region_FINAL.out   -common_args Manual_common_args.txt -native_pdb  mutate_G7-C8_A_site_1t0e_RNA_B.pdb -recreate_silent_struct True -skip_clustering True -extract_pdb False


##############################################################################################################################
START_argv=copy.deepcopy(argv)

no_graphic_string= parse_options( argv, "no_graphic", "" )

if(use_new_src_code()): 
	rna_swa_test_exe= get_rosetta_EXE_specified_no_graphic_string("swa_rna_main", no_graphic_string) 
else:
	rna_swa_test_exe = get_rosetta_EXE_specified_no_graphic_string("rna_swa_test", no_graphic_string)

database_folder= get_rosetta_database_folder() 

##############################################################################################################################


silent_files = parse_options( argv, "silent_file", [""] )

if(silent_files==[""]): error_exit_with_message("silent_files==[\"\"]")

silent_files_string=""
silent_files_string_underscore=""


for silent_file in silent_files:
	if(exists( silent_file )==False): error_exit_with_message("silent_file %s doesn't exist!" %(abspath(silent_file) ) )
	silent_files_string+=silent_file + " "
	if(silent_files_string_underscore==""):
		silent_files_string_underscore=silent_file
	else:
		silent_files_string_underscore+="_" + silent_file 


 
#remove_variant_types = parse_options( argv, "remove_variant_types", "True" )

cluster_rmsd= parse_options( argv, "cluster_rmsd", 2.0 )

###########################################################################################################################
user_distinguish_pucker = parse_options( argv, "distinguish_pucker", "" ) 

distinguish_pucker = ""

if(user_distinguish_pucker!=""):

	if(user_distinguish_pucker!="false" and user_distinguish_pucker!="true"): error_exit_with_message("Invalid user_distinguish_pucker (%s)" %(user_distinguish_pucker) )

	distinguish_pucker=user_distinguish_pucker

	print "user passed in distinguish_pucker option, distinguish_pucker=%s" %(distinguish_pucker)

else:

	if(cluster_rmsd>1.99999):
		print "cluster_rmsd (%s) > 1.999999 Setting distinguish_pucker to false!" %(cluster_rmsd)
		distinguish_pucker="false"
	else:
		print "cluster_rmsd (%s) < 1.999999 Setting distinguish_pucker to true!" %(cluster_rmsd)
		distinguish_pucker="true"

if(distinguish_pucker!="false" and distinguish_pucker!="true"): error_exit_with_message("Invalid distinguish_pucker (%s)" %(distinguish_pucker) )
##############################################################################################################################


add_lead_zero_to_tag=parse_options( argv, "add_lead_zero_to_tag", "true")

clusterer_rename_tags=parse_options( argv, "clusterer_rename_tags", "true")

skip_clustering= parse_options( argv, "skip_clustering", "False")

native_pdb= parse_options(argv, "native_pdb", "")
if(native_pdb!="" and exists(native_pdb)==False): error_exit_with_message("native_pdb %s doesn't exist!" %(native_pdb) )

extract_pdb= parse_options( argv, "extract_pdb", "True")

write_score_only= parse_options( argv, "write_score_only", "False")

align_only_over_base_atoms= parse_options( argv, "align_only_over_base_atoms", "True")

simple_append_map= parse_options( argv, "simple_append_map", "false") 

recreate_silent_struct=parse_options( argv, "recreate_silent_struct", "False" )

use_best_neighboring_shift_RMSD=parse_options( argv, "use_best_neighboring_shift_RMSD", "False" )

suite_cluster_OPTION=parse_options(argv, "suite_cluster_OPTION", "False")

suite_cluster_rmsd=parse_options(argv, "suite_cluster_rmsd", 0.00000)

if(suite_cluster_rmsd<0.00000): error_exit_with_message("suite_cluster_rmsd=(%s)<0.00000" %(suite_cluster_rmsd) )

if(suite_cluster_rmsd>0.00001): suite_cluster_OPTION=True	

if( (suite_cluster_OPTION==True) and (suite_cluster_rmsd<0.00001) ):
	suite_cluster_rmsd=cluster_rmsd

print "suite_cluster_OPTION=%s, suite_cluster_rmsd=%s " %(suite_cluster_OPTION,suite_cluster_rmsd)


num_pose_kept= parse_options( argv, "num_pose_kept", 1000 )

quick_alignment=parse_options(argv, "quick_alignment", "False") #Uncomment this option on Nov 5, 2011!

output_filename= parse_options( argv, "output_filename", "")

output_foldername = parse_options( argv, "output_foldername", "")

perform_VDW_rep_screen= parse_options(argv, "perform_VDW_rep_screen", "False")

perform_filters= parse_options(argv, "perform_filters", "False")

VDW_rep_screen_physical_pose_clash_dist_cutoff= parse_options(argv, "VDW_rep_screen_physical_pose_clash_dist_cutoff", "1.2")  #Was "0.12" before March 28, 2011...MISTAKE SHOULD BE "1.2"

if(VDW_rep_screen_physical_pose_clash_dist_cutoff[0]=="N"): VDW_rep_screen_physical_pose_clash_dist_cutoff="-" + VDW_rep_screen_physical_pose_clash_dist_cutoff[1:]

full_length_loop_rmsd_clustering= parse_options( argv, "full_length_loop_rmsd_clustering", "False")

ignore_FARFAR_no_auto_bulge_tag= parse_options( argv, "ignore_FARFAR_no_auto_bulge_tag", "False") #For post processing!

ignore_FARFAR_no_auto_bulge_parent_tag= parse_options( argv, "ignore_FARFAR_no_auto_bulge_parent_tag", "False") #For post processing!

ignore_unmatched_virtual_res= parse_options( argv, "ignore_unmatched_virtual_res", "False") #For post processing!

redirect_out_log= parse_options( argv, "redirect_out_log", "True") #For post processing!

VERBOSE= parse_options(argv, "VERBOSE", "false")


common_args=""
common_args_file = parse_options(argv, "common_args", "")
if(common_args_file != ""):
	print "common_args_file= ", common_args_file
	common_args=open( common_args_file  ).readlines()[0].strip()

	
print common_args

###After parsing there should only be one argv left which is the python file name: SWA_rna_build_dagman.py####
if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )


print
print


keep_pose_in_memory=True #Dec 11, 2011

if(skip_clustering): keep_pose_in_memory=False #Dec 11, 2011

if(use_best_neighboring_shift_RMSD): keep_pose_in_memory=True #Dec 11, 2011



#Optional



cluster_rmsd_int=int(cluster_rmsd*100)

suite_cluster_rmsd_int=int(suite_cluster_rmsd*100)

if(output_filename==""):

	if(use_best_neighboring_shift_RMSD==True): 

		output_filename="neighbor_shift_RMSD_%d_" %(cluster_rmsd_int)

		if(suite_cluster_OPTION==True): output_filename+="suite_%s_" %(suite_cluster_rmsd_int)

		output_filename+="%s" %(silent_files_string_underscore.replace('/','_') )

	else:
		output_filename='cluster_%d_%s' % (cluster_rmsd_int , silent_files_string_underscore.replace('/','_') )


		if(recreate_silent_struct==True): output_filename="recreate_" + output_filename

		if(suite_cluster_OPTION==True): output_filename="suite_clustering_" + output_filename

		if(perform_filters):
			output_filename='ROSETTA_filter_%s' %(silent_files_string_underscore.replace('/','_') )


		if(skip_clustering and recreate_silent_struct):
			output_filename='recal_rmsd_%s' %(silent_files_string_underscore.replace('/','_') )

log_file=append_to_basename('CLUSTER_LOG_' , output_filename)


print "cluster_rmsd=%f Angstrom" % cluster_rmsd
print "output_filename= %s" %(output_filename)

if(exists(output_filename)):
	print "output_filename (%s) already exist! .... removing... " %(output_filename) 
	submit_subprocess( 'rm %s'  % output_filename )

###########################################################################
command=rna_swa_test_exe
command += " -algorithm rna_cluster"

if(common_args==""):
	command += " -whole_struct_cluster_radius %s " % cluster_rmsd
else:
	if(suite_cluster_OPTION==True):
		command += " -suite_cluster_radius %s " %(suite_cluster_rmsd)

	command += " -loop_cluster_radius %s " %(cluster_rmsd)
	command += " %s " %(common_args)


command += " -in:file:silent " + silent_files_string

if(option_name_exist(command, "in:file:silent_struct_type")):

	old_args=get_option_name_args_safe(command, "in:file:silent_struct_type", "")

	if(old_args!="binary_rna"): error_exit_with_message('old_args!="binary_rna"')

else:
	command += " -in:file:silent_struct_type  binary_rna"

command += " -database %s " %(database_folder)
command += " -out:file:silent " + output_filename
command += " -add_lead_zero_to_tag %s " %(add_lead_zero_to_tag)
command += " -clusterer_two_stage_clustering false " 
command += " -clusterer_rename_tags %s " %(clusterer_rename_tags)
command += " -simple_append_map %s " %(simple_append_map)
command += " -VERBOSE %s " %(VERBOSE)

if(align_only_over_base_atoms==False): command += " -clusterer_align_only_over_base_atoms false "

if(recreate_silent_struct==True): command += " -recreate_silent_struct true "

if(use_best_neighboring_shift_RMSD==True): command += " -clusterer_use_best_neighboring_shift_RMSD true "

if(distinguish_pucker=="false"): command += " -distinguish_pucker false "

if(skip_clustering==True): command += " -skip_clustering true "

command += " -clusterer_num_pose_kept %d " %(num_pose_kept)

if(native_pdb!=""): command += ' -native %s ' %(native_pdb)

if(write_score_only==True): command += ' -clusterer_write_score_only true '

if(quick_alignment==True): command += ' -clusterer_quick_alignment true '

if(perform_VDW_rep_screen==True): 
	command += '-clusterer_perform_VDW_rep_screen true '
	command += '-VDW_rep_screen_physical_pose_clash_dist_cutoff %s ' %(VDW_rep_screen_physical_pose_clash_dist_cutoff)

if(perform_filters==True):
	command += '-clusterer_perform_filters true '

if(full_length_loop_rmsd_clustering==True):
	command += "-clusterer_full_length_loop_rmsd_clustering true " 

if(ignore_FARFAR_no_auto_bulge_tag==True): command += "-clusterer_ignore_FARFAR_no_auto_bulge_tag true "

if(ignore_FARFAR_no_auto_bulge_parent_tag==True): command += "-clusterer_ignore_FARFAR_no_auto_bulge_parent_tag true "

if(ignore_unmatched_virtual_res==True): command += "-clusterer_ignore_unmatched_virtual_res true "

if(keep_pose_in_memory==False): command += "-clusterer_keep_pose_in_memory false " #Dec 11, 2011

if(redirect_out_log): command += ' > %s' %(log_file)

print command 

ensure_no_duplicate_options(command)

submit_subprocess( command )

if(output_foldername==""):
	output_foldername="pose_" + output_filename.replace('/','_')

if(exists(output_foldername)==True): submit_subprocess("rm -r %s " %(output_foldername))


if(extract_pdb):
	extract_command='SWA_extract_pdb.py -silent_file %s -move_silent_file_into_folder True -output_foldername %s' % (output_filename, output_foldername)
	if(no_graphic_string!=""): extract_command+=' -no_graphic %s ' %(no_graphic_string)

	print extract_command
	submit_subprocess( extract_command )
	submit_subprocess(  'mv %s %s/.' % (output_filename,output_foldername) )

	if(redirect_out_log): submit_subprocess(  'mv %s %s/.' % (log_file, output_foldername) )


print "-------------------------------------------------------------------------------------------"
print "Successfully RAN: %s" %(list_to_string(START_argv))
print "-------------------------------------------------------------------------------------------"



