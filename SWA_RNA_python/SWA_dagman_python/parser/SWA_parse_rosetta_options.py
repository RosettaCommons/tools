#!/usr/bin/env python

######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_parse_options import parse_options, parse_seq_num_list_option
######################################################################


def get_general_rosetta_common_args_option(argv): #These are shared with FARFAR!

	common_args = ' '

	#ROSETTA TRUNK (Dec 2011):
	#option.add( basic::options::OptionKeys::in::file::silent_read_through_errors, "will ignore decoys with errors and continue reading" ).def(false);
	#option.add( basic::options::OptionKeys::out::file::output_virtual, "Output virtual atoms in output of PDB" ).def(false);
	#option.add( basic::options::OptionKeys::in::file::silent_struct_type, "Type of SilentStruct object to use in silent-file input" ).def("protein");

	common_args += ' -in:file:silent_struct_type binary_rna '

	if(is_release_mode()==False):	#Consistency check!

		common_args += ' -output_virtual true '
		common_args += ' -silent_read_through_errors false '  #silent_read_through_errors used to be true for clusterer! #move this option to common_args on Sept 21, 2011
		common_args += ' -output_extra_RMSDs true ' #For extra check.

	force_field_file=parse_options( argv, "force_field_file", "stepwise/rna/rna_hires_07232011_with_intra_base_phosphate.wts" )
	force_syn_chi_res_list= parse_options( argv, "force_syn_chi_res_list", [ -1 ] ) #May 01, 2011
	force_north_ribose_list=parse_options(argv, "force_north_ribose_list", [-1] ) #May 03, 2011
	force_south_ribose_list=parse_options(argv, "force_south_ribose_list", [-1] ) #May 03, 2011
	protonated_H1_adenosine_list=parse_options(argv, "protonated_H1_adenosine_list", [-1] ) #May 03, 2011
	rna_torsion_potential_folder=parse_options(argv, "rna_torsion_potential_folder", "" ) #May 06, 2011

	native_virtual_res = parse_options( argv, "native_virtual_res", [ -1 ] )


	if(force_field_file==""): error_exit_with_message("User need to specify force_field_file!")

	common_args += ' -score:weights %s ' %(force_field_file) #April 9th, 2011

	if(len(force_syn_chi_res_list)>0): common_args += ' -force_syn_chi_res_list %s ' %( list_to_string(force_syn_chi_res_list) )

	if(len(force_north_ribose_list)>0): common_args += ' -force_north_ribose_list %s ' %( list_to_string(force_north_ribose_list) )

	if(len(force_south_ribose_list)>0): common_args += ' -force_south_ribose_list %s ' %( list_to_string(force_south_ribose_list) )

	if(len(protonated_H1_adenosine_list)>0): common_args += ' -protonated_H1_adenosine_list %s ' %( list_to_string(protonated_H1_adenosine_list) )

	if(rna_torsion_potential_folder!=""):
		if(use_new_src_code()):
			common_args += ' -score:rna_torsion_potential %s ' %(rna_torsion_potential_folder)
		else:
			common_args += ' -score:rna_torsion_folder %s ' %(rna_torsion_potential_folder)

	if(len(native_virtual_res) > 0):	common_args += ' -native_virtual_res %s ' %(list_to_string(native_virtual_res) )

	return common_args

################################################################################

def get_rosetta_common_args_option(argv):

	common_args = ' '

	common_args += get_general_rosetta_common_args_option(argv)

	FAST=parse_options( argv, "fast", "false" )
	PARIN_FAVORITE_OUTPUT = parse_options( argv, "parin_favorite_output", "true" )

	fixed_res = parse_seq_num_list_option( argv, "fixed_res" )
	minimize_res = parse_options( argv, "minimize_res", [ -1 ] )
	terminal_res = parse_options( argv, "terminal_res", [ -1 ] )
	virtual_res = parse_options( argv, "virtual_res", [ -1 ]  )
	rmsd_res = parse_seq_num_list_option( argv, "rmsd_res" )
	native_alignment_res = parse_seq_num_list_option( argv, "native_alignment_res" )
	jump_point_pair_list = parse_options( argv, "jump_point_pair_list", [""] )
	alignment_res_list = parse_options( argv, "alignment_res_list", [""]  )
	fa_stack_base_base_only=parse_options( argv, "fa_stack_base_base_only", "true" )
	allow_chain_boundary_jump_partner_right_at_fixed_BP= parse_options( argv, "allow_chain_boundary_jump_partner_right_at_fixed_BP", "false"  )  #Hacky...need this to get square RNA to work! Nov 6, 2010
	allow_fixed_res_at_moving_res= parse_options( argv, "allow_fixed_res_at_moving_res", "false"  )  #Hacky...need this to get Hermann Duplex working!, Nov 15, 2010

	#Move VDW rep option from samplerer_args to common_args on March 20, 2011. This is for post_process VDW_rep_screening in SWA_clusterer
	VDW_rep_screen_info = parse_options( argv, "VDW_rep_screen_info", [""]  )
	VDW_rep_alignment_RMSD_CUTOFF = parse_options( argv, "VDW_rep_alignment_RMSD_CUTOFF", "0.001"  ) #Nov 12, 2010, to allow for imperfect alignment
	VDW_rep_delete_matching_res = parse_options( argv, "VDW_rep_delete_matching_res", [""]  ) #delete residues in VDW_rep_pose that exist in the working_pose
	VDW_rep_optimize_memory_usage = parse_options( argv, "VDW_rep_optimize_memory_usage", "false")


	if(FAST=="true"): common_args += ' -fast true'

	if (PARIN_FAVORITE_OUTPUT=="false"): error_exit_with_message('PARIN_FAVORITE_OUTPUT=="false"!')

	if( len( fixed_res )==0 and len( minimize_res )==0):
		error_exit_with_message("User need to specify either fixed_res or minimize_res")

	if(len( fixed_res ) > 0): common_args += ' -fixed_res %s ' %(list_to_string_with_dashes(fixed_res) )

	if(len(minimize_res) > 0): common_args += ' -minimize_res %s ' %(list_to_string_with_dashes(minimize_res) )

	if(len(terminal_res) > 0): common_args += ' -terminal_res %s ' %(list_to_string_with_dashes(terminal_res) )

	if(len(virtual_res) > 0): common_args += ' -virtual_res %s ' %(list_to_string_with_dashes(virtual_res) )

	if(len( rmsd_res ) > 0): common_args += ' -rmsd_res %s ' %(list_to_string_with_dashes(rmsd_res) )

	if(jump_point_pair_list!=[""]):

		if(use_new_src_code()):
			common_args += ' -jump_point_pairs %s ' %(list_to_string(jump_point_pair_list))
		else:
			common_args += ' -fixed_BP_list %s ' %(list_to_string(jump_point_pair_list))

	if(alignment_res_list!=[""]): common_args += ' -alignment_res %s ' %(list_to_string(alignment_res_list))

	if( len(native_alignment_res)>0 ):
			common_args += ' -native_alignment_res %s ' %(list_to_string_with_dashes(native_alignment_res))

	if(fa_stack_base_base_only=="false"): common_args += ' -score::fa_stack_base_base_only false '

	if(allow_chain_boundary_jump_partner_right_at_fixed_BP=="true"): common_args += ' -allow_chain_boundary_jump_partner_right_at_fixed_BP true '

	if(allow_fixed_res_at_moving_res=="true"): common_args += ' -allow_fixed_res_at_moving_res true '

	########################################################################################
	if(VDW_rep_screen_info!=[""]):
		check_valid_VDW_rep_screen_info_list(VDW_rep_screen_info)

		VDW_rep_screen_info_command= ' -VDW_rep_screen_info '
		for k in VDW_rep_screen_info: VDW_rep_screen_info_command += '%s ' % k
		print VDW_rep_screen_info_command
		common_args += VDW_rep_screen_info_command

		if(VDW_rep_alignment_RMSD_CUTOFF!="0.001"):
			common_args += ' -VDW_rep_alignment_RMSD_CUTOFF %s ' %(VDW_rep_alignment_RMSD_CUTOFF)

		if(VDW_rep_delete_matching_res!=[""]):
			common_args += ' -VDW_rep_delete_matching_res %s ' %(list_to_string(VDW_rep_delete_matching_res) )

	if(VDW_rep_optimize_memory_usage=="true"): common_args += ' -VDW_rep_optimize_memory_usage true '
	########################################################################################

	return common_args

##############################################################################################################

def get_rosetta_samplerer_args_option(argv):

	native_rmsd_screen=parse_options( argv, "native_rmsd_screen", "false" )
	rmsd_screen=parse_options( argv, "rmsd_screen", 0.0 )
	sampler_native_screen_rmsd_cutoff=parse_options( argv, "sampler_native_screen_rmsd_cutoff", "2.0" )
	sampler_num_pose_kept = parse_options( argv, "sampler_num_pose_kept", 108 )
	sampler_extra_epsilon_rotamer = parse_options( argv, "sampler_extra_epsilon_rotamer", "true" ) #change default to true on April 17th, 2011
	sampler_extra_beta_rotamer = parse_options( argv, "sampler_extra_beta_rotamer", "false" )
	sampler_extra_anti_chi_rotamer = parse_options( argv, "sampler_extra_anti_chi_rotamer", "false" ) #split to anti and syn on Juen 17th, 2011
	sampler_extra_syn_chi_rotamer = parse_options( argv, "sampler_extra_syn_chi_rotamer", "false" )   #split to anti and syn on Juen 17th, 2011
	sampler_include_torsion_value_in_tag = parse_options( argv, "sampler_include_torsion_value_in_tag", "true" )
	allow_base_pair_only_centroid_screen = parse_options( argv, "allow_base_pair_only_centroid_screen", "false" ) #this only effect the dinucleotide floating base mode..#change to false on April 9th, 2011
	do_not_sample_multiple_virtual_sugar = parse_options( argv, "do_not_sample_multiple_virtual_sugar", "false" ) #this only effect the dinucleotide floating base mode..
	sample_ONLY_multiple_virtual_sugar = parse_options( argv, "sample_ONLY_multiple_virtual_sugar", "false" ) #this only effect the dinucleotide floating base mode..
	allow_bulge_at_chainbreak = parse_options( argv, "allow_bulge_at_chainbreak", "true"  ) #combine_long_file filterer need to know this information (still need to implement)??? #change to true on April 9th, 2011

	tether_jump = parse_options( argv, "tether_jump", "true")

	include_syn_chi = parse_options( argv, "include_syn_chi", "true"  ) #Move from common_args on Oct 28, 2011

	sampler_perform_o2star_pack = parse_options( argv, "sampler_perform_o2star_pack", "true" ) #Move from common_args on Oct 28, 2011

	sampling_args=' '

	if (native_rmsd_screen=="true"):
		sampling_args += ' -sampler_native_rmsd_screen true '
		sampling_args += ' -sampler_native_screen_rmsd_cutoff %s ' %(sampler_native_screen_rmsd_cutoff) #move this into if loop on April 9th, 2011

	if( rmsd_screen > 0.0 ): sampling_args += ' -rmsd_screen %d ' % (rmsd_screen)

	if(sampler_num_pose_kept!=108): sampling_args += ' -sampler_num_pose_kept %d ' %(sampler_num_pose_kept)

	if(sampler_extra_epsilon_rotamer=="false"): sampling_args += ' -sampler_extra_epsilon_rotamer false '

	if(sampler_extra_beta_rotamer=="true"): sampling_args += ' -sampler_extra_beta_rotamer true '

	if(sampler_extra_anti_chi_rotamer=="true"): sampling_args += ' -sampler_extra_anti_chi_rotamer true '

	if(sampler_extra_syn_chi_rotamer=="true"): sampling_args += ' -sampler_extra_syn_chi_rotamer true '

	if(sampler_include_torsion_value_in_tag=="false"): sampling_args += ' -sampler_include_torsion_value_in_tag false '

	if(allow_base_pair_only_centroid_screen=="true"): sampling_args += ' -allow_base_pair_only_centroid_screen true '

	if(do_not_sample_multiple_virtual_sugar=="true"): sampling_args += ' -do_not_sample_multiple_virtual_sugar true '

	if(sample_ONLY_multiple_virtual_sugar=="true"): sampling_args += ' -sample_ONLY_multiple_virtual_sugar true '

	if(allow_bulge_at_chainbreak=="false"): sampling_args += ' -allow_bulge_at_chainbreak false '

	if(tether_jump=="false"): sampling_args += ' -tether_jump false '

	if(include_syn_chi=="false"): sampling_args += ' -include_syn_chi false'

	if(sampler_perform_o2star_pack=="false"): sampling_args +=  ' -sampler_perform_o2star_pack false'

	return sampling_args

##############################################################################################################

def get_rosetta_clusterer_args_option(argv):

	suite_cluster_radius = parse_options( argv, "suite_cluster_radius", "1.0"  )
	loop_cluster_radius = parse_options( argv, "loop_cluster_radius", "0.7"  )
	clusterer_quick_alignment=parse_options( argv, "clusterer_quick_alignment", "false" )
	clusterer_optimize_memory_usage=parse_options( argv, "clusterer_optimize_memory_usage", "false" )
	clusterer_keep_pose_in_memory=parse_options( argv, "clusterer_keep_pose_in_memory", "true" )
	clusterer_two_stage_clustering=parse_options( argv, "clusterer_two_stage_clustering", "false" ) #Change default to true on Oct 11, 2010 #change to false on April 9th, 2011


	cluster_args=' '

	cluster_args += ' -suite_cluster_radius %s ' % suite_cluster_radius
	cluster_args += ' -loop_cluster_radius %s ' % loop_cluster_radius

	if(clusterer_quick_alignment=="true"):
		cluster_args += ' -clusterer_quick_alignment true '

	if(clusterer_optimize_memory_usage=="true"):
		cluster_args += ' -clusterer_optimize_memory_usage true '

	if(clusterer_keep_pose_in_memory=="false"):
		cluster_args += ' -clusterer_keep_pose_in_memory false '

	if(clusterer_two_stage_clustering=="true"):
		cluster_args += ' -clusterer_two_stage_clustering true '

	return cluster_args

