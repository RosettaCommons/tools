#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
import os
import copy
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################

##########################10:40 PM Feb 06 2012: UUAC_tetraloop###########
#SWA_copy_non_full_length_regions.py -num_elements 5

#rsync -avz   /Volumes/Sept_2011_Ext_HD_Parin/minirosetta/11_to_12_2011_CS_BLIND/Fox_UNAC_loop/Nov_10_SWA_UNAC_tetraloop/main_folder/UUAC_tetraloop/PARTIAL_REGIONS_FOLDER/ .

##########################10:40 PM Feb 06 2012: UGAC_tetraloop###########
#SWA_copy_non_full_length_regions.py -num_elements 5

#rsync -avz   /Volumes/Sept_2011_Ext_HD_Parin/minirosetta/11_to_12_2011_CS_BLIND/Fox_UNAC_loop/Nov_10_SWA_UNAC_tetraloop/main_folder/UGAC_tetraloop/PARTIAL_REGIONS_FOLDER/ .

##########################10:37 PM Feb 06 2012: UAAC_tetraloop###########
#SWA_copy_non_full_length_regions.py -num_elements 5

#rsync -avz   /Volumes/Sept_2011_Ext_HD_Parin/minirosetta/11_to_12_2011_CS_BLIND/Fox_UNAC_loop/Nov_10_SWA_UNAC_tetraloop/main_folder/UAAC_tetraloop/PARTIAL_REGIONS_FOLDER/ .

##########################10:38 PM Feb 06 2012: UCAC_tetraloop###########
#SWA_copy_non_full_length_regions.py -num_elements 5


#rsync -avz   /Volumes/Sept_2011_Ext_HD_Parin/minirosetta/11_to_12_2011_CS_BLIND/Fox_UNAC_loop/Nov_10_SWA_UNAC_tetraloop/main_folder/UCAC_tetraloop/PARTIAL_REGIONS_FOLDER/ .

##########################10:33 PM Feb 06 2012: TRNA_YE_MET##############
#SWA_copy_non_full_length_regions.py -num_elements 8

#rsync -avz   /Volumes/Sept_2011_Ext_HD_Parin/minirosetta/11_to_12_2011_CS_BLIND/Nikonowicz_TRNA/Dec_20_SWA_TRNA_ANTI_CODON_PART_2/main_folder/TRNA_YE_MET/PARTIAL_REGIONS_FOLDER/ .

##########################10:30 PM Feb 06 2012 : UUAAGU_hexaloop##############
#SWA_copy_non_full_length_regions.py -num_elements 7

#rsync -avz   /Volumes/Sept_2011_Ext_HD_Parin/minirosetta/10_2011/10_2011_UUAAGU/Oct_22_SWA_UUAAGU_hexaloop_ALL_PATHS/main_folder/UUAAGU_Hexaloop/PARTIAL_REGIONS_FOLDER/ .

##########################Feb 05 2012 : CUUG Tetraloop##############
#SWA_copy_non_full_length_regions.py -num_elements 5

#rsync -avz   /Volumes/Sept_2011_Ext_HD_Parin/minirosetta/11_2011/Nov_03_SWA_CUUG_loop/main_folder/CUUG_tetraloop/PARTIAL_REGIONS_FOLDER/ .

###########################Feb_09 2012: 1LNT#######################################
#SWA_copy_non_full_length_regions.py -num_elements 10 -cutpoint_open 5 -missing_allow_bulge_right_next_to_input_helix True -remove_old_RMSD_cols True > LOG_SWA_copy_non_full_length_regions.txt


#rsync -avz   /Volumes/Sept_2011_Ext_HD_Parin/minirosetta/10_2011/Oct_20_1LNT_SWA_RUN/main_folder/SRP_HELIX_8/PARTIAL_REGIONS_FOLDER/ .

##########################09:36 AM Feb 09 2012 : UUAAGU_hexaloop##############
#SWA_copy_non_full_length_regions.py -num_elements 7 -remove_old_RMSD_cols True  > LOG_SWA_copy_non_full_length_regions.txt

#rsync -avz   /Volumes/Sept_2011_Ext_HD_Parin/minirosetta/10_2011/10_2011_UUAAGU/Oct_22_SWA_UUAAGU_hexaloop_ALL_PATHS/main_folder/UUAAGU_Hexaloop/PARTIAL_REGIONS_FOLDER/ .


##########################06:46 AM Feb 10 2012 : GGAC_duplex##############
#SWA_copy_non_full_length_regions.py -num_elements 6 -cutpoint_open 3 -missing_allow_bulge_right_next_to_input_helix True  > LOG_SWA_copy_non_full_length_regions.txt

#rsync -avz   /Volumes/Sept_2011_Ext_HD_Parin/minirosetta/11_2011/Nov_02_SWA_GGAC_duplex/main_folder/GGAC/PARTIAL_REGIONS_FOLDER/ .

##########################06:46 AM Feb 10 2012 : GGAC_duplex##############
#SWA_copy_non_full_length_regions.py -num_elements 7 -cutpoint_open 3 -missing_allow_bulge_right_next_to_input_helix False  > LOG_SWA_copy_non_full_length_regions.txt


#rsync -avz   /Volumes/Sept_2011_Ext_HD_Parin/minirosetta/11_to_12_2011_CS_BLIND/Sigel_39mer/Nov_27_SWA_Sigel_39mer/main_folder/Sigel_39mer/PARTIAL_REGIONS_FOLDER/ .


##########################08:01 PM Feb 14 2012 : GGAC_duplex##############
#SWA_copy_non_full_length_regions.py -num_elements 5 -cutpoint_open 2 -missing_allow_bulge_right_next_to_input_helix False  > LOG_SWA_copy_non_full_length_regions.txt


#rsync -avz  /Volumes/Sept_2011_Ext_HD_Parin/minirosetta/11_to_12_2011_CS_BLIND/Schawble_HARF1/Nov_12_SWA_HARF1/main_folder/short_Chimp_HAR1/PARTIAL_REGIONS_FOLDER/ .
###################################################################

copy_argv=copy.deepcopy(argv)

#########################################################
missing_allow_bulge_right_next_to_input_helix=parse_options( argv, "missing_allow_bulge_right_next_to_input_helix", "False" )

remove_old_RMSD_cols=parse_options( argv, "remove_old_RMSD_cols", "False" )

num_elements=parse_options( argv, "num_elements", 0 )

if(num_elements<=0): error_exit_with_message("num_elements=(%s)==0" %(num_elements))

cutpoint_open=parse_options( argv, "cutpoint_open", 0 ) #if cutpoint_open!=0, then double-stranded internal loop! #THIS REFERS TO THE ELEMENT BEFORE THE CUTPOINT OPEN

PARTIAL_REGIONS_FOLDER="PARTIAL_REGIONS_FOLDER/" 

if(exists(PARTIAL_REGIONS_FOLDER)): error_exit_with_message("PARTIAL_REGIONS_FOLDER/ (%s) already exist!" %(PARTIAL_REGIONS_FOLDER))

submit_subprocess("mkdir %s" %(PARTIAL_REGIONS_FOLDER))

all_region_pair_list=[]

if(cutpoint_open==0): #Single-stranded hairpin loop

	for lower_element in range(0, num_elements): #if num_elements=4 -> 0, 1, 2, 3

		upper_element_bound=lower_element

		if(upper_element_bound==0): upper_element_bound=num_elements

		for upper_element in range(0, upper_element_bound):

			if(lower_element==0 and upper_element==0): continue

			all_region_pair_list.append([lower_element, upper_element])
			

else:

	#Lower-left square/Build inward
	for lower_element_ID in range(cutpoint_open+1, num_elements+1): #if num_elements=10, cutpoint_open=5 -> 6, 7, 8, 9, 0

		lower_element=lower_element_ID

		if(lower_element==num_elements): lower_element=0
		
		for upper_element in range(0, cutpoint_open):  #if num_elements=10, cutpoint_open=5 -> 0, 1, 2, 3, 4

			if(lower_element==0 and upper_element==0): continue

			all_region_pair_list.append([lower_element, upper_element])	


	#Upper-right square/Build outward
	for lower_element in range(1, cutpoint_open+1): #if num_elements=10, cutpoint_open=5 -> 1, 2, 3, 4, 5

		if(lower_element==num_elements): lower_element=0
		
		for upper_element in range(cutpoint_open, num_elements): #if num_elements=10, cutpoint_open=5 -> 5, 6, 7, 8, 9

			if(lower_element==cutpoint_open and upper_element==cutpoint_open): continue

			all_region_pair_list.append([lower_element, upper_element])	

	if(missing_allow_bulge_right_next_to_input_helix==False):
		
		################Lower-left square/Build inward##################
		if( (cutpoint_open+2) > num_elements):
			print "----------------------------------\n----------------------------------"
			print "(cutpoint_open+2) > num_elements!"
			print "----------------------------------\n----------------------------------"
		elif( (cutpoint_open+2) == num_elements):
			all_region_pair_list.append([0, cutpoint_open])

		else:
			all_region_pair_list.append([cutpoint_open+2, cutpoint_open])
		#################################################################

		if( (cutpoint_open-2) < 0):
			print "----------------------------------\n----------------------------------"
			print "cutpoint_open-2) < 0!"
			print "----------------------------------\n----------------------------------"
		else:
			all_region_pair_list.append([cutpoint_open, cutpoint_open-2])
		#################################################################

		################Upper-right square/Build outward#################
		if( 2 > cutpoint_open):
			print "----------------------------------\n----------------------------------"
			print "2 > cutpoint_open!"
			print "----------------------------------\n----------------------------------"
		else:
			all_region_pair_list.append([2, 0])

		#################################################################
		if( (num_elements-2) < cutpoint_open):
			print "----------------------------------\n----------------------------------"
			print "(num_elements-2) < cutpoint_open!"
			print "----------------------------------\n----------------------------------"
		else:
			all_region_pair_list.append([0, num_elements-2])
		#################################################################



################################################################################################

partial_region_silent_file_list=[]

for region_pair in all_region_pair_list:

	if(len(region_pair)!=2): error_exit_with_message("len(region_pair)!=2")

	lower_element=region_pair[0]
	upper_element=region_pair[1]

	silent_file="region_%s_%s_sample.cluster.out" %(lower_element, upper_element)

	if(exists(silent_file)==False): error_exit_with_message("silent_file (%s) doesn't exist!" %(silent_file))

	####Screen out full regions####
	Is_full_length_region=False

	if(lower_element==0):
		if(upper_element==(num_elements-1)): Is_full_length_region=True
	else:
		if(upper_element==(lower_element-1)): Is_full_length_region=True

	print "region_pair=%s |silent_file=%s --->Is_full_length_region=%s" %(list_to_string(region_pair), silent_file, Is_full_length_region)

	if(cutpoint_open!=0 and Is_full_length_region): error_exit_with_message("cutpoint_open!=0 BUT Is_full_length_region==TRUE!")

	if(Is_full_length_region==False): partial_region_silent_file_list.append(silent_file)

################################################################################################

partial_region_silent_file_list.sort()

for silent_file in partial_region_silent_file_list:

	copy_command="cp %s %s/%s" %(silent_file, PARTIAL_REGIONS_FOLDER, silent_file)

	print "copy_command= %s" %(copy_command)

	submit_subprocess(copy_command)




if(remove_old_RMSD_cols):

	print "BEFORE ENTER PARTIAL_REGIONS_FOLDER | os.path.abspath(\".\")=",  os.path.abspath(".")

	print "os.chdir( PARTIAL_REGIONS_FOLDER )"
	os.chdir( PARTIAL_REGIONS_FOLDER )

	print "AFTER ENTER PARTIAL_REGIONS_FOLDER | os.path.abspath(\".\")=",   os.path.abspath(".")

	for silent_file in partial_region_silent_file_list:

		if(exists(silent_file)==False): error_exit_with_message("silent_file (%s) doesn't exist!" %(silent_file))	

		replace_scoreline_silent_file="replace_scoreline_%s" %(silent_file)

		if(exists(replace_scoreline_silent_file)==True): error_exit_with_message("replace_scoreline_silent_file (%s) already exist!" %(replace_scoreline_silent_file))	

		replace_silent_scoreline_command="replace_silent_scoreline.py -keep_column_names  score     fa_atr     fa_rep    fa_intra_rep    fa_intra_RNA_base_phos_atr    fa_intra_RNA_base_phos_rep    lk_nonpolar    lk_nonpolar_intra_RNA    hack_elec_rna_phos_phos    ch_bond    rna_torsion    rna_sugar_close    fa_stack    hbond_sr_bb_sc    hbond_lr_bb_sc    hbond_sc    hbond_intra    geom_sol_intra_RNA    CI_geom_sol    atom_pair_constraint    angle_constraint    rna_bulge    linear_chainbreak    all_rms       O_rmsd    O_loop_rmsd    O_V_rms    O_V_loop_rms    O_PBP_rmsd  -infile %s " %(silent_file) 

		submit_subprocess(replace_silent_scoreline_command)

		if(exists(replace_scoreline_silent_file)==False): error_exit_with_message("replace_scoreline_silent_file (%s) doesn't exist!" %(replace_scoreline_silent_file))	

		remove_command="rm %s " %(silent_file)

		move_command="mv %s %s" %(replace_scoreline_silent_file, silent_file)

		print "remove_command=%s | move_command=%s" %(remove_command, move_command)

		submit_subprocess(remove_command)

		submit_subprocess(move_command)

print "----------------------------------------------------------------------------------------------------------------------------"
print "Sucessfully RAN: %s " %(copy_argv)
print "----------------------------------------------------------------------------------------------------------------------------"

