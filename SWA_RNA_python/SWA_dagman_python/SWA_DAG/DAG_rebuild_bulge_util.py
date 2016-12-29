#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser,abspath
from sys import exit, argv
import string
from time import sleep
import glob
import time
import fileinput
import os
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options, replace_arg_value
######################################################################

######################################################################
HOMEDIR = expanduser('~') 

 

if(use_new_src_code()): 
	RNA_SWA_MAIN_EXE  = get_rosetta_EXE("swa_rna_main") 
	RNA_SWA_UTIL_EXE  = get_rosetta_EXE("swa_rna_util") 
else:
	RNA_SWA_MAIN_EXE  = get_rosetta_EXE("rna_swa_test") 
	RNA_SWA_UTIL_EXE  = get_rosetta_EXE("parin_test") 

DATABASE_FOLDER= get_rosetta_database_folder()

print "RNA_SWA_MAIN_EXE =%s" %(RNA_SWA_MAIN_EXE)
print "RNA_SWA_UTIL_EXE =%s" %(RNA_SWA_UTIL_EXE)
print "DATABASE_FOLDER=%s" %(DATABASE_FOLDER)

##################################################################################################################
def parse_sampling_options(argv, native_pdb):

	fine_sampling = parse_options( argv, "fine_sampling", "True" ) 

	freeze_chi_o2star = parse_options( argv, "freeze_chi_o2star", "False" )

	sampling_options=""

	if(fine_sampling):
		sampling_options+= ' -sampler_cluster_rmsd 0.5 '  #ensure that we keep a diverse range of conformations
	else:
		sampling_options+= ' -sampler_cluster_rmsd 3.0 '  #ensure that we keep a diverse range of conformations

	if(freeze_chi_o2star==True): 
		error_exit_with_message("Jan 24, 2012: freeze_chi_o2star mode have not yet been tested!")
		#Jan 24, 2012: STILL NEED TO TEST THIS MODE (especially if using the REP only force-field since IF NOT the 2'-OH hydrogen will be mis-positioned)!!
		sampling_options+= " -sampler_perform_o2star_pack false "
		sampling_options+= " -minimizer_perform_o2star_pack false " 


	if(native_pdb!=""): sampling_options+=" -native " + native_pdb

	return sampling_options

################################################################################################################################################################
def extract_sampling_options_from_common_args(common_args_file, virtual_res_list, location, start_folder):

	sampling_options=" "

	common_args=safe_readlines( common_args_file )[0].split()

	####BUGGY ####include_syn_chi = "False" #Oct 22, 2011: Notice that "False" is not a correct Rosetta c++ option. But at least on local mac laptop this leads to the desire outcome in that syn_chi is NOT sampled! CONFIRMED THIS BEHAVOIR ON Jan 23, 2012 on local machine as well.

	include_syn_chi="false" #Change to this on Jan 23, 2012!

	force_field_file = parse_options(common_args, "score:weights", "")

	#Note that if use force_field->06_23_2011_REP_only_rna_loop_hires.wts" then need to set freeze_chi_o2star to true!

	protonated_H1_adenosine_list= parse_options( common_args, "protonated_H1_adenosine_list", [0] )
	force_syn_chi_res_list= parse_options( common_args, "force_syn_chi_res_list", [-1] )
	force_north_ribose_list= parse_options( common_args, "force_north_ribose_list", [-1] )
	force_south_ribose_list= parse_options( common_args, "force_south_ribose_list", [-1] )
	VDW_rep_screen_info=parse_options( common_args , "VDW_rep_screen_info", [""])
	VDW_rep_delete_matching_res=parse_options( common_args , "VDW_rep_delete_matching_res", "false")

	#TRIED ADDING THIS on Nov 20, 2011 | BECAREFUL! THIS ARG IS INDEPENTENTLY PARSED OUTSIDE in the main DAG_rebuild_bulge.py script.
	#HOWEVER DECIDE TO COMMENT OUT SINCE PREFER THE SIMPLE FOLD_TREE WHEN REBUILDING THE BULGE. SIMPLE FOLD TREE SHOULD BE MORE ROBUST IN HANDLING THE COMPLEX MOTIF OF FARFAR CASES
	#cutpoint_open_list= parse_options( common_args, "cutpoint_open", [0]) 


	global_sample_res_list=parse_options( common_args, "global_sample_res_list", [0]) #BECAREFUL! THIS ARG IS INDEPENTLY PARSED OUTSIDE in the main DAG_rebuild_bulge.py script.
	#####################################################################################################################

	force_syn_chi_res_list=list_intersection(force_syn_chi_res_list, virtual_res_list)
	force_north_ribose_list=list_intersection(force_north_ribose_list, virtual_res_list)
	force_south_ribose_list=list_intersection(force_south_ribose_list, virtual_res_list)

	#####################################################################################################################


	if(force_field_file!=""): sampling_options+=  " -score:weights %s " %(force_field_file)

	sampling_options+=  " -include_syn_chi %s " %(include_syn_chi) 


	if(len(protonated_H1_adenosine_list)>0): sampling_options+= " -protonated_H1_adenosine_list %s " %(list_to_string(protonated_H1_adenosine_list)) 

	if(len(force_syn_chi_res_list)>0):			 sampling_options+= ' -force_syn_chi_res_list %s ' %(list_to_string(force_syn_chi_res_list) ) 

	if(len(force_north_ribose_list)>0): 		 sampling_options+= ' -force_north_ribose_list %s ' %(list_to_string(force_north_ribose_list) )

	if(len(force_south_ribose_list)>0): 		 sampling_options+= ' -force_south_ribose_list %s ' %(list_to_string(force_south_ribose_list) )

	if(VDW_rep_screen_info!=[""]):

		current_folder=os.path.abspath(".")
		os.chdir( start_folder )
		
		for n in range(len(VDW_rep_screen_info)):
			if(VDW_rep_screen_info[n][-4:]!=".pdb"): continue

			if(location=="LOCAL"): VDW_rep_screen_info[n]="../" + VDW_rep_screen_info[n]
	
			VDW_rep_screen_info[n]= os.path.abspath(VDW_rep_screen_info[n])

		sampling_options+= " -VDW_rep_screen_info %s " %(list_to_string(VDW_rep_screen_info))
		sampling_options+= " -VDW_rep_delete_matching_res %s " %(VDW_rep_delete_matching_res)

		os.chdir( current_folder )

	#if(len(cutpoint_open_list)>0): sampling_options+= " -cutpoint_open %s " %(list_to_string(cutpoint_open_list))

	if(len(global_sample_res_list)==0): error_exit_with_message("len(global_sample_res_list)==0")
	sampling_options+= " -global_sample_res_list %s " %(list_to_string(global_sample_res_list))

	return sampling_options

################################################################################################################################################################
def run_rebuild_bulge_sampling(input_silent_file, input_tag, rebuild_seq_num, cutpoint_closed, user_fixed_res_list, total_res, fasta_file, location, sampling_options):

	#Feb 18, 2012: NOTE, missing cutpoint_open (list)!!!

	tag_ID="RES_%d_CC_%d" %(rebuild_seq_num, cutpoint_closed)

	output_silent_file="REBUILD_BULGE_RES_%s_CC_%s.out" %(rebuild_seq_num , cutpoint_closed)

	sampling_argv= RNA_SWA_MAIN_EXE + ' -algorithm rna_sample ' 

	sampling_argv+= " -database %s " %(DATABASE_FOLDER)

	sampling_argv+= " -in:file:silent_struct_type binary_rna " #Feb 18, 2012.

	sampling_argv+= " -output_virtual true "

	sampling_argv+= " -silent_read_through_errors false " #Feb 18, 2012.

	sampling_argv+= " -rebuild_bulge_mode true " #IMPORTANT

	sampling_argv+= " -fasta %s " %(fasta_file)

	sampling_argv+= " -input_res %s " %( list_to_string(range(1, total_res+1)) )

	#sampling_argv+= " -fa_stack_base_base_only false "######COMMENT OUT ON OCT 23, 2011: Should not have this line since fa_stack_base_base_only should be true!!!

	sampling_argv+= " -PBP_clustering_at_chain_closure true "
	 
	sampling_argv+= " -reinitialize_CCD_torsions false "
	
	sampling_argv+= " -allow_bulge_at_chainbreak true "

	sampling_argv+= " -sampler_extra_epsilon_rotamer true " #Nov 30, 2010

	sampling_argv+= " -sampler_extra_beta_rotamer false " #Nov 30, 2010

	sampling_argv+= " -sampler_extra_syn_chi_rotamer false " #Nov 30, 2010

	sampling_argv+= " -sampler_extra_anti_chi_rotamer false " #Nov 30, 2010

	if(location=="LOCAL"):
		sampling_argv+= " -VERBOSE true"
		sampling_argv+= " -output_pdb true "
	else:
		sampling_argv+= " -VERBOSE false"
		sampling_argv+= " -output_pdb false "
	
	sampling_argv+= " -alignment_res 1-%d " %(total_res)
	sampling_argv+= " -fixed_BP_list 1-%d " %(total_res)

	#########################Specific_Argv##################################
	specific_argv=""

	fixed_res_list=[]

	if(len(user_fixed_res_list)>0):

		fixed_res_list=user_fixed_res_list
		print "Use USER_input fixed_res_list: %s " %(fixed_res_list)	

	else:

		for fix_seq_num in range(1,total_res+1):	 

			if(fix_seq_num==rebuild_seq_num): continue

			fixed_res_list.append(fix_seq_num)

		print "Use default fixed_res_list= %s " %(list_to_string(fixed_res_list))

	specific_argv+= " -fixed_res %s " %(list_to_string(fixed_res_list))

	specific_argv+= " -rmsd_res %d " %(rebuild_seq_num)

	specific_argv+= " -tags %s " %(input_tag) 

	specific_argv+= " -in:file:silent %s " %(input_silent_file)

	specific_argv+= " -out:file:silent %s " %(output_silent_file) 

	#specific_argv+= " -global_sample_res_list %d " %(rebuild_seq_num) ##Comment out on Oct 24, 2011! See parse_sampling_options() for replacement!

	specific_argv+= " -sample_res %d " %(rebuild_seq_num)

	specific_argv+= " -cutpoint_closed %d " %(cutpoint_closed)

	command = sampling_argv + ' ' + specific_argv + ' '  + sampling_options + ' > rosetta_sampling_outfile_%s.txt 2> rosetta_sampling_errfile_%s.txt' %(tag_ID, tag_ID)

	print  '\n', "command= %s" %(command), '\n'

	sys.stdout.flush()
	sys.stderr.flush()

	if(location=="LOCAL"): #graphic sometime clash the code.
		system(command)
	else:
		submit_subprocess( command ) 

	sys.stdout.flush()
	sys.stderr.flush()

	return output_silent_file
	##########################################################
	

################################################################################
def extract_start_tag_info(start_silent_file, start_tag):

	start_desc_column=get_description_column(start_silent_file)

	found_start_score_line=False


	start_score_line=""
	start_column_line=""
	ANNOTATED_SEQUENCE=""

	num_start_score_line=0
	num_start_column_line=0

	for line in open(start_silent_file):

		if(found_start_score_line==False):

			if(line.count('SCORE:') == 0): continue #Line contain the word SCORE

			if(line[0:6]!="SCORE:"): error_exit_with_message("line.count('SCORE:') != 0 but line[0:6]!=\"SCORE:\"") 
			
			if(line.count('description') != 0): 
			
				start_column_line=line
				num_start_column_line+=1
				continue

			else:

				if( (line.split()[start_desc_column])==start_tag ): 
					start_score_line=line
					found_start_score_line=True

					num_start_score_line+=1
		
		else:
			if(line.count('ANNOTATED_SEQUENCE') != 0 ): 
				ANNOTATED_SEQUENCE=line
				found_start_score_line=False


	if(num_start_column_line!=1): error_exit_with_message("num_start_column_line=(%s)!=1 for start_silent_file (%s) !" %(num_start_column_line,start_silent_file))

	if(num_start_score_line!=1): error_exit_with_message("num_start_score_line=(%s)!=1 for start_tag (%s) !" %(num_start_score_line, start_tag))

	if(ANNOTATED_SEQUENCE==""): error_exit_with_message("ANNOTATED_SEQUENCE=\"\"")

	if(start_score_line==""): error_exit_with_message("start_score_line==\"\"")

	if(start_column_line==""): error_exit_with_message("start_column_line==\"\"")

	if( (start_score_line.split()[start_desc_column])!=start_tag ): error_exit_with_message("(start_score_line.split()[start_desc_column])!=start_tag")

	if( (start_column_line.split()[start_desc_column])!="description" ): error_exit_with_message("(start_column_line.split()[start_desc_column])!=\"description\"")


	ANNOTATED_SEQUENCE=ANNOTATED_SEQUENCE.split()[1]

	return (ANNOTATED_SEQUENCE, start_score_line, start_column_line)


################################################################################################################################################################


def parse_annotated_sequence(ANNOTATED_SEQUENCE):

	seq_list=[]
	virtual_res_list=[]
	cutpoint_lower_list=[]
	cutpoint_upper_list=[]

	seq_num=1
	n=0

	#ANNOTATED_SEQUENCE: g[RGU_p:Virtual_Phosphate]ccga[RAD_p:Virtual_RNA_Residue]a[RAD_p:rna_cutpoint_lower_p:Virtual_Phosphate_p:Virtual_RNA_Residue_Upper]a[RAD_p:rna_cutpoint_upper]ggc S_7
	while(n<len(ANNOTATED_SEQUENCE)):

		char=ANNOTATED_SEQUENCE[n]

		if(char=='['): 
			end_bracket=ANNOTATED_SEQUENCE.find(']', n)	
	
			Is_virtual_res=False

			if( ANNOTATED_SEQUENCE.count('Virtual_RNA_Residue:', n, end_bracket)>0  ): Is_virtual_res=True

			if( ANNOTATED_SEQUENCE.count('Virtual_RNA_Residue]', n, end_bracket+1)>0 ): Is_virtual_res=True

			if( ANNOTATED_SEQUENCE.count('Virtual_RNA_Residue_p:', n, end_bracket)>0  ): Is_virtual_res=True

			if( ANNOTATED_SEQUENCE.count('Virtual_RNA_Residue_p]', n, end_bracket+1)>0 ): Is_virtual_res=True

			Is_cutpoint_lower=(ANNOTATED_SEQUENCE.count('rna_cutpoint_lower', n, end_bracket+1)>0 ) #something weird rna_cutpoint_lower_p (sometime have _p, sometime don't)
			Is_cutpoint_upper=(ANNOTATED_SEQUENCE.count('rna_cutpoint_upper', n, end_bracket+1)>0 )

			if(Is_virtual_res): virtual_res_list.append(seq_num-1)
			if(Is_cutpoint_lower): cutpoint_lower_list.append(seq_num-1)
			if(Is_cutpoint_upper): cutpoint_upper_list.append(seq_num-1)

			n=end_bracket+1
		else:
			seq_list.append(char )
			seq_num=seq_num+1
			n=n+1

	#Check 
	if(len(cutpoint_lower_list)!=len(cutpoint_upper_list)): error_exit_with_message("len(cutpoint_lower_list)!=len(cutpoint_upper_list)" )
	for cutpoint_ID in range(len(cutpoint_lower_list)):
		if((cutpoint_lower_list[cutpoint_ID]+1)!=cutpoint_upper_list[cutpoint_ID]):
			error_exit_with_message("cutpoint_lower_list[%s]+1)=%s=!=%s=cutpoint_upper_list[%s]" %(cutpoint_ID, cutpoint_lower_list[cutpoint_ID], cutpoint_upper_list[cutpoint_ID],cutpoint_ID) )

	return (seq_list, virtual_res_list, cutpoint_lower_list)


################################################################################################################################################################
def figure_out_cutpoint_closed_list(main_algorithm, rebuild_seq_num_list, start_silent_file, start_tag, total_res):

	if(len(rebuild_seq_num_list)==0): #Early return
		cutpoint_closed_list=[]
		print "-----------------------------------------------"
		print "(len(rebuild_seq_num_list)==0) ...EARLY RETURN!"
		print "rebuild_seq_num_list= ", rebuild_seq_num_list
		print "cutpoint_closed_list= ", cutpoint_closed_list
		print "-----------------------------------------------"
		return cutpoint_closed_list


	cutpoint_closed_file="rebuild_bulge_cutpoint_closed_list.txt"

	if(exists(cutpoint_closed_file)): error_exit_with_message("cutpoint_closed_file (%s) already exist!" %(cutpoint_closed_file))

	command= RNA_SWA_UTIL_EXE + ' -algorithm figure_out_rebuild_bulges_cutpoint_closed_list ' 

	command+= " -database %s " %(DATABASE_FOLDER)

	command+= " -virtual_res %s " %(list_to_string(rebuild_seq_num_list))
	
	command+= " -tags %s " %(start_tag) 

	command+= " -in:file:silent %s " %(start_silent_file)

	command+= " -output_text_file %s " %(cutpoint_closed_file)
	
	command+= " -base_sampler_algorithm %s " %(main_algorithm)

	submit_subprocess(command)

	if(exists(cutpoint_closed_file)==False): error_exit_with_message("cutpoint_closed_file (%s) doesn't exist!" %(cutpoint_closed_file))

	line_list=open(cutpoint_closed_file).readlines()

	if(len(line_list)!=1): error_exit_with_message("len(line_list)=(%s)!=1" %(line_list))

	cutpoint_closed_list=line_list[0].split()

	if(cutpoint_closed_list[0]!="-cutpoint_closed_list"): 
		error_exit_with_message("cutpoint_closed_list[0]!=\"-cutpoint_closed_list\", for cutpoint_closed_line (%s)" %(line_list[0]))

	cutpoint_closed_list=cutpoint_closed_list[1:]

	cutpoint_closed_list=	map( lambda x : int(x) , cutpoint_closed_list )

	if(len(cutpoint_closed_list)!=len(rebuild_seq_num_list)): 
		print "ERROR rebuild_seq_num_list=%s" %(list_to_string(rebuild_seq_num_list))
		print "ERROR cutpoint_closed_list=%s" %(list_to_string(cutpoint_closed_list))
		print "ERROR line_list[0]=%s" %(line_list[0])
		error_exit_with_message("len(cutpoint_closed_list)!=len(rebuild_seq_num_list)")

	for bulge_ID in range(len(rebuild_seq_num_list)):

		valid_cutpoint_pos=False

		cutpoint_closed=cutpoint_closed_list[bulge_ID]
		virtual_res=rebuild_seq_num_list[bulge_ID] 

		if(cutpoint_closed not in [virtual_res-1, virtual_res]):
			print "ERROR bulge_ID=%s" %(bulge_ID)
			print "ERROR rebuild_seq_num_list=%s" %(list_to_string(rebuild_seq_num_list))
			print "ERROR cutpoint_closed_list=%s" %(list_to_string(cutpoint_closed_list))
			print "ERROR line_list[0]=%s" %(line_list[0])
			error_exit_with_message("cutpoint_closed (%s) not in [virtual_res-1 (%s), virtual_res (%s)]" %(cutpoint_closed, virtual_res-1, virtual_res) )

		if((cutpoint_closed+1)>total_res):
			error_exit_with_message("(cutpoint_closed+1)(%s)>total_res(%s)", (cutpoint_closed+1,total_res) )


	print "rebuild_seq_num_list= ", rebuild_seq_num_list
	print "cutpoint_closed_list= ", cutpoint_closed_list

	return cutpoint_closed_list

################################################################################################################################################################

'''
def figure_out_cutpoint_closed_list_old(main_algorithm, rebuild_seq_num_list, cutpoint_open_list, cutpoint_lower_list, total_res):

#########START OF AUTO_REBUILD_RES_AND_CUTPOINT==True################: 
	cutpoint_closed_list=[] #clear the list

	for virtual_res in rebuild_seq_num_list

		##determine the corresponding cutpoint_closed location.
		print "Figuring out corresponding_cutpoint_closed to rebuild_bulge_res= %d " %(virtual_res) 


		if(main_algorithm=="FARFAR"):
			print 'main_algorithm=="FARFAR"'
			cutpoint_closed_list.append(virtual_res)
		else: #SWA
			if(len(cutpoint_open_list)<=1):
				if(len(cutpoint_lower_list)!=1):
					error_exit_with_message("len(cutpoint_open_list)==1 but len(cutpoint_lower_list)!=1, cutpoint_lower_list=%s" %(list_to_string(cutpoint_lower_list) ) )

			cutpoint_lower=cutpoint_lower_list[0]

			#Virtual could be right at the cutpoint
			if(virtual_res in cutpoint_lower_list):
				print "CASE: virtual_res in cutpoint_lower_list"
				cutpoint_closed_list.append(virtual_res)
			elif((virtual_res-1) in cutpoint_lower_list):
				print "CASE: (virtual_res-1) in cutpoint_lower_list"
				cutpoint_closed_list.append(virtual_res-1)
			else: 
				print "DINCULEOTIDE/FLOATING BASE CASES:"

				if(len(cutpoint_open_list)>1): #multiple strand like VCII junction..

					print "CASE: len(cutpoint_open_list)>1 , guessing cutpoint_closed location."
					cutpoint_closed_list.append(virtual_res) #not enough info to figure out actual cutpoint_closed location

				elif(len(cutpoint_open_list)==1):

					print "CASE: len(cutpoint_open_list)== 1"
					cutpoint_open=cutpoint_open_list[0]

					num_closing_helix_BP=2

					close_helix_location="" #either close chain at upper or lower element.

					####Nov 25, 2011:Case #3 and Case #4 are for the allow_bulge_right_next_to_input_helix mode.
					####Nov 25, 2011: Note that there are situation where the condition below will incorrectly choose the close_helix_location.
					####For example if one strand is empty then cutpoint_lower then Case #1 and Case #2 are true.
					####The problem is more prevalent when allow_bulge_right_next_to_input_helix mode=True:
					####When one strand contain 1 nt, now more than one of Case #1-4 can be true
					####Even when one strand contain 2 nts, now both Case #3 and Case #4 can both be true.
					####The most general solution is to check (inside Rosetta), which of the two possible cutpoint_close does not have the idealize geometry and close there!

					if( ( num_closing_helix_BP in cutpoint_lower_list) or ( (total_res-num_closing_helix_BP) in cutpoint_lower_list) ): #Case 1
						close_helix_location="LOWER"
					elif( ( (cutpoint_open-num_closing_helix_BP) in cutpoint_lower_list) or ( (cutpoint_open+num_closing_helix_BP) in cutpoint_lower_list) ): #Case2
						close_helix_location="UPPER"
					elif( ( (num_closing_helix_BP+1) in cutpoint_lower_list) or ( (total_res-num_closing_helix_BP-1) in cutpoint_lower_list) ): #Case 3 (Nov 25, 2011)
						close_helix_location="LOWER"
					elif( ( (cutpoint_open-num_closing_helix_BP-1) in cutpoint_lower_list) or ( (cutpoint_open+num_closing_helix_BP+1) in cutpoint_lower_list) ): #Case #4 (Nov 25, 2011)
						close_helix_location="UPPER"
					else:
						error_exit_with_message("could not figure out close_helix_location" )

					Is_prepend=False						
					if(close_helix_location=="LOWER"):
						if(virtual_res<= cutpoint_open):
							Is_prepend=True
						else:
							Is_prepend=False
					elif(close_helix_location=="UPPER"):
						if(virtual_res<= cutpoint_open):				
							Is_prepend=False
						else:
							Is_prepend=True
					else:
						error_exit_with_message("SHOULD NOT REACH THIS POINT IN THE CODE!" )
				
					if(Is_prepend):
						cutpoint_closed_list.append(virtual_res-1)	
					else:
						cutpoint_closed_list.append(virtual_res)	

				else: #loop region.

					print "CASE: len(cutpoint_open_list)== 0"

					if(virtual_res>=cutpoint_lower): #prepend
						cutpoint_closed_list.append(virtual_res-1)	
					else:
						cutpoint_closed_list.append(virtual_res)	

	#########END OF AUTO_REBUILD_RES_AND_CUTPOINT==True################: 


	print "rebuild_seq_num_list= ", rebuild_seq_num_list
	print "cutpoint_closed_list= ", cutpoint_closed_list

	if(len(rebuild_seq_num_list)!=len(cutpoint_closed_list)): error_exit_with_message("len(rebuild_seq_num_list)!=len(cutpoint_closed_list)" )

	for bulge_ID in range(len(rebuild_seq_num_list)):
		cutpoint_closed=cutpoint_closed_list[bulge_ID]
		rebuild_seq_num=rebuild_seq_num_list[bulge_ID]

		if((cutpoint_closed+1)>total_res):
			error_exit_with_message("(cutpoint_closed+1)(%s)>total_res(%s)", (cutpoint_closed+1,total_res) )

		if(cutpoint_closed!=rebuild_seq_num and cutpoint_closed!=(rebuild_seq_num-1)):
			 error_exit_with_message("cutpoint_closed!=rebuild_seq_num and cutpoint_closed!=(rebuild_seq_num-1)")

	return cutpoint_closed_list
'''
################################################################################

def sort_silent_file_by_score(input_silent_file):

	sorted_silent_file="sorted_%s" %(basename(input_silent_file))

	if(exists(sorted_silent_file)): 
		error_exit_with_message("sorted_silent_file (%s) already exist!" %(sorted_silent_file) )

	cluster_command="SWA_cluster.py -silent_file %s -output_filename %s -cluster_rmsd 0.0 -extract_pdb False -redirect_out_log False " %(input_silent_file, sorted_silent_file)

	#print "------------------------------------------------------------------------------------------------------------------"
	#print "cluster_command= %s" %(cluster_command) 
	#print "------------------------------------------------------------------------------------------------------------------"

	submit_subprocess( cluster_command )

	return sorted_silent_file
