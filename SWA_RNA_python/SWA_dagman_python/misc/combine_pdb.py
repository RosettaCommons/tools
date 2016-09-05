#!/usr/bin/env python

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
from os.path import basename, dirname, exists, expanduser
from time import sleep
import os
from os import popen 
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################

remove_storage_folder=True

main_folder=os.path.abspath(".")


#combine_pdb.py aligned_HIV_TAR_upper_helix_to_DAS_FARFAR_01.pdb-1-2 DAS_FARFAR_01.pdb-1-10 aligned_HIV_TAR_upper_helix_to_DAS_FARFAR_01.pdb-7-8 -output_pdb_name HIV_TAR_upper_element_include_G11_C25.pdb
   

#combine_pdb.py 2aw7_RNA_A.pdb-913-919 2aw7_RNA_A.pdb-1389-1392 2aw7_RNA_A.pdb-11-14  -output_pdb_name   mock_native_excited_A_site_2aw7.pdb

#combine_pdb.py rosetta_3nj6_RNA_A.pdb-1-10 rosetta_3nj6_partner_RNA_A.pdb-1-10 -output_pdb_name 3nj6_duplex.pdb


#combine_pdb.py 3p59_RNA_A.pdb-98-98 3p59_RNA_A.pdb-80-85 3p59_RNA_A.pdb-61-62 3p59_RNA_A.pdb-86-88 3p59_RNA_A.pdb-93-97 -output_pdb_name corner_GH_3p59_RNA_A.pdb


#combine_pdb.py 3p59_RNA_A.pdb-25-25 3p59_RNA_A.pdb-5-15 3p59_RNA_A.pdb-20-24 -output_pdb_name corner_AB_3p59_RNA_A.pdb


#combine_pdb.py 3p59_RNA_A.pdb-50-50 3p59_RNA_A.pdb-30-40 3p59_RNA_A.pdb-45-49 -output_pdb_name corner_CD_3p59_RNA_A.pdb


#combine_pdb.py 3p59_RNA_A.pdb-75-75 3p59_RNA_A.pdb-55-65 3p59_RNA_A.pdb-70-74 -output_pdb_name corner_EF_3p59_RNA_A.pdb

#COMBINE HIV TAR APICAL LOOP AND UPPER IDEALIZE HELIX.

#combine_pdb.py DASLAB_TS002_1.pdb-86-100 DASLAB_TS002_1.pdb-1-10 DASLAB_TS002_1.pdb-71-85 DASLAB_TS002_1.pdb-11-20 DASLAB_TS002_1.pdb-56-70 DASLAB_TS002_1.pdb-21-30 DASLAB_TS002_1.pdb-41-50 DASLAB_TS002_1.pdb-53-55 DASLAB_TS002_1.pdb-31-40 -strict_line_length_check False -output_pdb_name RENUM_DASLAB_TS002_1.pdb


#combine_pdb.py DASLAB_TS002_2.pdb-86-100 DASLAB_TS002_2.pdb-1-10 DASLAB_TS002_2.pdb-71-85 DASLAB_TS002_2.pdb-11-20 DASLAB_TS002_2.pdb-56-70 DASLAB_TS002_2.pdb-21-30 DASLAB_TS002_2.pdb-41-50 DASLAB_TS002_2.pdb-53-55 DASLAB_TS002_2.pdb-31-40 -strict_line_length_check False -output_pdb_name RENUM_DASLAB_TS002_2.pdb


#combine_pdb.py DASLAB_TS002_3.pdb-86-100 DASLAB_TS002_3.pdb-1-10 DASLAB_TS002_3.pdb-71-85 DASLAB_TS002_3.pdb-11-20 DASLAB_TS002_3.pdb-56-70 DASLAB_TS002_3.pdb-21-30 DASLAB_TS002_2.pdb-41-50 DASLAB_TS002_3.pdb-53-55 DASLAB_TS002_3.pdb-31-40 -strict_line_length_check False -output_pdb_name RENUM_DASLAB_TS002_3.pdb

#combine_pdb.py DASLAB_TS002_4.pdb-86-100 DASLAB_TS002_4.pdb-1-10 DASLAB_TS002_4.pdb-71-85 DASLAB_TS002_4.pdb-11-20 DASLAB_TS002_4.pdb-56-70 DASLAB_TS002_4.pdb-21-30 DASLAB_TS002_4.pdb-41-50 DASLAB_TS002_4.pdb-53-55 DASLAB_TS002_4.pdb-31-40 -strict_line_length_check False -output_pdb_name RENUM_DASLAB_TS002_4.pdb


#combine_pdb.py aligned_HIV_TAR_upper_helix_to_DAS_FARFAR_03.pdb-1-2 DAS_FARFAR_03.pdb-1-10 aligned_HIV_TAR_upper_helix_to_DAS_FARFAR_03.pdb-7-8 -output_pdb_name HIV_TAR_upper_element_include_G11_C25.pdb
   
#combine_pdb.py DASLAB_TS002_5.pdb-86-100 DASLAB_TS002_5.pdb-1-10 DASLAB_TS002_5.pdb-71-85 DASLAB_TS002_5.pdb-11-20 DASLAB_TS002_5.pdb-56-70 DASLAB_TS002_5.pdb-21-30 DASLAB_TS002_5.pdb-41-50 DASLAB_TS002_5.pdb-53-55 DASLAB_TS002_5.pdb-31-40 -strict_line_length_check False -output_pdb_name RENUM_DASLAB_TS002_5.pdb



#2PN4
#combine_pdb.py res_3_CC_2_S_000004.pdb-1-4 res_6_CC_5_S_000013.pdb-5-16 -output_pdb_name S_0_res_3_CC_2_S_000004_res_6_CC_5_S_000013.pdb
      
#Hermann phase 2
#combine_pdb.py 283d_Single_strand.pdb-1-12 283d_sym_partner_RNA_A.pdb-1-12  -output_pdb_name 283d_duplex.pdb

#VC_II 4 nucleotides loop
#combine_pdb.py RES_6_CC_6_S_000000.pdb-1-7  RES_8_CC_8_S_000000.pdb-8-10 -output_pdb_name Rebuild_bulge_FARFAR_CAUA_loop_S_0.pdb


#INVERT
#combine_pdb.py 2PN3_Square_RNA.pdb-12-17  2PN3_Square_RNA.pdb-1-11 -output_pdb_name consistency_check_invert_2PN3_Square_RNA.pdb

#INVERT
#combine_pdb.py REGION_9_0_S_000000.pdb-1-3  REGION_0_1_S_000000.pdb-4-7  REGION_9_0_S_000000.pdb-9-9 -output_pdb_name SWA_minimized_idealized_BP_START_INVERT_2PN3_Square_RNA.pdb                       


#combine_pdb.py WITH_MOCK_TRAIL_2_S_000002.pdb-1-1 WITH_MOCK_TRAIL_2_S_000002.pdb-11-17 -output_pdb_name TRAIL_2_S_000002.pdb

#combine_pdb.py S_000003.pdb-1-1  ROTATE_MINUS_ONE_2PN3_Square_RNA.pdb-2-10 S_000003.pdb-2-2 S_000003.pdb-4-9 -output_pdb_name delete_G17_S_000003_mock_res_2_10.pdb
#combine_pdb.py S_000005.pdb-1-1  ROTATE_MINUS_ONE_2PN3_Square_RNA.pdb-2-10 S_000005.pdb-2-2 S_000005.pdb-4-9 -output_pdb_name delete_G12_S_000005_mock_res_2_10.pdb


#MINUS_ONE
#combine_pdb.py S_000005.pdb-17-17  S_000005.pdb-1-16 -output_pdb_name ROTATE_MINUS_ONE_Nov_6_S_000005.pdb

#combine_pdb.py 2PN3_Square_RNA.pdb-17-17  2PN3_Square_RNA.pdb-1-16 -output_pdb_name ROTATE_MINUS_ONE_2PN3_Square_RNA.pdb

#combine_pdb.py REGION_0_2_S_000000.pdb-1-2 REGION_0_2_S_000000.pdb-4-8 -output_pdb_name REGION_0_1_S_000000.pdb

#combine_pdb.py REGION_10_0_S_000000.pdb-1-1 REGION_10_0_S_000000.pdb-3-8 -output_pdb_name REGION_11_0_S_000000.pdb



#combine_pdb.py REGION_0_1_S_000000.pdb-1-4  REGION_9_0_S_000000.pdb-6-9 -output_pdb_name SWA_minimized_idealized_BP_START_rotate_seq_2PN3_Square_RNA.pdb-1-3                        

#combine_pdb.py 2PN3_Square_RNA.pdb-15-17 2PN3_Square_RNA.pdb-1-11 2PN3_Square_RNA.pdb-12-14

#combine_pdb.py inner_strand_2PN3_Square_RNA.pdb-4-6 inner_strand_2PN3_Square_RNA.pdb-1-3 

#combine_pdb.py FIXED_AU_BP_aligned_to_AB_corner.pdb-1-1 FIXED_CG_BP_aligned_to_CG_corner.pdb-2-2 short_corner_AB_1_10.pdb-1-6
#combine_pdb.py FIXED_A4U5_aligned_upper_VDW_rep_screener.pdb-1-8 FIXED_C4G5_aligned_lower_VDW_rep_screener.pdb-1-8   #note that I mistakenly switched lower and upper VDW_screen_info from the START.

#combine_pdb.py S_0_seq_num_7_cc_7_S_000009.pdb-1-11 S_0_seq_num_17_cc_16_S_000015.pdb-12-22



user_input_output_pdb_name= parse_options( argv, "output_pdb_name", "")

ABS_PATH_user_input_output_pdb_name=os.path.abspath(user_input_output_pdb_name)

strict_line_length_check =parse_options( argv, "strict_line_length_check", "True")

consistency_check =parse_options( argv, "consistency_check", "False")

combine_pdb_folder=parse_options( argv, "combine_pdb_folder", "")

if(combine_pdb_folder==""):
	combine_pdb_folder="storage_combine_pdb/"

if(exists(combine_pdb_folder)): 
	print "%s already exist...removing..." %(combine_pdb_folder)
	submit_subprocess("rm -r %s " %(combine_pdb_folder) )
 
submit_subprocess( "mkdir %s" %(combine_pdb_folder) )



pdb_segment_info_list=[]

output_pdb="combine"

for n in range(1,len(argv)):
	print argv[n]
	pdb_segment_split=argv[n].split('-')

	output_pdb+="_%s" %(argv[n].replace("-","_"))

	pdb_file=pdb_segment_split[0]

	pdb_segment_split=pdb_segment_split[1:]
	
	if(exists(pdb_file)==False): error_exit_with_message("pdb_file (%s) doesn't exist" %(pdb_file) ) 

	copy_pdb_file="%s/%s" %(combine_pdb_folder, basename(pdb_file) )
	submit_subprocess( "cp %s %s" %(pdb_file, copy_pdb_file) )
	pdb_file=basename(pdb_file)


	pdb_segment_info={}

	if((len(pdb_segment_split) % 2) != 0): error_exit_with_message("len(pdb_segment_split) % 2 != 0. pdb_segment_split=%s" %(list_to_string(pdb_segment_split)) )



	for n in range(len(pdb_segment_split)):
		try:
			pdb_segment_split[n]=int(pdb_segment_split[n])
		except:
			error_exit_with_message("cannot convert pdb_segment_split[%d] to int, pdb_segment_split[%d]=%s" %(n,n, pdb_segment_split[n]) )

	for n in range( 1, len(pdb_segment_split) ):

		if( pdb_segment_split[n]<pdb_segment_split[n-1] ): 
			error_exit_with_message("pdb_segment_split[%d]<pdb_segment_split[%d], pdb_segment_split=%s" %(n, n-1, list_to_string(pdb_segment_split) ) )

	
	pdb_segment_info["pdb_file"]=pdb_file
	pdb_segment_info["segment_list"]=pdb_segment_split

	pdb_segment_info_list.append(pdb_segment_info)

output_pdb+=".pdb"



os.chdir( combine_pdb_folder )



for n in range(len(pdb_segment_info_list)):
	print "pdb_segment %d: " %(n+1) , pdb_segment_info_list[n] 
	if(exists(pdb_segment_info_list[n]["pdb_file"])==False): error_exit_with_message("pdb_file (%s) doesn't exist" %(pdb_segment_info_list[n]["pdb_file"]) ) 

if(user_input_output_pdb_name!=""):
	output_pdb=ABS_PATH_user_input_output_pdb_name

print "output_pdb= %s" %(output_pdb)

if(exists(output_pdb)):

	if(user_input_output_pdb_name!=""):
		error_exit_with_message("user_input_output_pdb_name (%s) already exist!" %(user_input_output_pdb_name))
	else:	
		print "Warning %s already exist...removing!" %(output_pdb)
		submit_subprocess( "rm %s " %(output_pdb) )
	

curr_res_num = 0
curr_atom_num = 0
last_resnum=""
resnum=""

OUPUT_PDB = open( output_pdb, 'w')

#ATOM      1  P    rA A   1       6.703 -48.708  19.088  1.00 21.73        
#0123456789012345678901234567890123456789012345678901234567890123456789
#0					1					2				 3					 4					5					6

for n in range(len(pdb_segment_info_list)):

	pdb_segment_info=pdb_segment_info_list[n]

	pdb_file=pdb_segment_info["pdb_file"]
	
	pdbslice_command="SWA_pdbslice.py %s -segments" %(pdb_file)
	sliced_pdb="sliced_%s_segments" %(pdb_file.replace(".pdb", "") )

	for seq_num in pdb_segment_info["segment_list"]:
		pdbslice_command+=" %d" %(seq_num)
		sliced_pdb+="_%d" %(seq_num)	

	sliced_pdb+=".pdb"
	pdbslice_command+=" %s" %(sliced_pdb)

	print "pdbslice_command=%s" %(pdbslice_command)

		
	submit_subprocess(pdbslice_command )
	if(exists(sliced_pdb)==False): error_exit_with_message("sliced_pdb %(s) doesn't exist!" %(slice_pdb))			 

	pdb_string_list=open(sliced_pdb).readlines()
	new_file=True

	for line in pdb_string_list:

		if(line=="END\n"): 
			#OUPUT_PDB.write("END\n")
			continue

		if(strict_line_length_check):
			if(len(line)!= 81): error_exit_with_message("len(line)!= 81, len(line)=%d, line=%s" %(len(line), line) )  
		else:
			if(len(line)< 80): error_exit_with_message("len(line)< 0, len(line)=%d, line=%s" %(len(line), line) )  


		if(line[0:4] != 'ATOM'): error_exit_with_message("line[0:4] != 'ATOM', line=%s" %(line) ) 

		if(line[16]!=' '): error_exit_with_message("line[16]!=' ' != 'ATOM', line=%s" %(line) ) 

		curr_atom_num += 1

		resnum = line[22:26]
		if((resnum!=last_resnum) or (new_file==True)):
			curr_res_num+=1

		new_file=False
		last_resnum = resnum

		line = '%s%5d%s%s%s' % (line[0:6],curr_atom_num,line[11:22], '%4d' % curr_res_num, line[26:] )
		OUPUT_PDB.write(line)
	######new on Nov 11, 2010#########
	#submit_subprocess("renumber_pdb_in_place.py %s" %(sliced_pdb) )
	##################################

	#submit_subprocess("cat %s >> %s " %(sliced_pdb, output_pdb))

OUPUT_PDB.close()
submit_subprocess("replace_chain_inplace.py %s A " %(output_pdb) )

#######################consistency check#######################################
if(consistency_check):
	renumber_pdb_in_place="renumber_pdb_in_place_" + basename(output_pdb)

	#if(dirname(output_pdb)!=""):
	#	renumber_pdb_in_place= os.path.abspath(dirname(output_pdb)) +'/' + renumber_pdb_in_place

	#print "renumber_pdb_in_place= %s" %(renumber_pdb_in_place)


	submit_subprocess("cp %s %s " %(output_pdb,renumber_pdb_in_place))
	submit_subprocess("renumber_pdb_in_place.py %s" %(renumber_pdb_in_place) )


	diff_line_list = popen('diff %s %s ' %(output_pdb, renumber_pdb_in_place) ).readlines()

	if(len(diff_line_list)!=0):
		for diff_line in diff_line_list:
			print "diff_line: %s" %(diff_line)
		error_exit_with_message("len(diff_line_list)!=0")


	#######################consistency check#######################################

	rosetta_ready_pdb=output_pdb.lower().replace(".pdb", "_RNA_A.pdb")

	if(output_pdb==rosetta_ready_pdb) : error_exit_with_message("output_pdb(%s)==rosetta_ready_pdb(%s)" %(output_pdb,rosetta_ready_pdb) )

	if(exists(rosetta_ready_pdb)==True): error_exit_with_message("rosetta_ready_pdb (%s) already exist!!" %(rosetta_ready_pdb))  

	submit_subprocess("SWA_make_rna_rosetta_ready.py %s" %(output_pdb))

	if(exists(rosetta_ready_pdb)==False): error_exit_with_message("rosetta_ready_pdb (%s) doesn't exist!!" %(rosetta_ready_pdb))  

	diff_line_list = popen('diff %s %s ' %(output_pdb, rosetta_ready_pdb) ).readlines()

	if(len(diff_line_list)!=0):
		for diff_line in diff_line_list:
			print "diff_line: %s" %(diff_line)
		error_exit_with_message("len(diff_line_list)!=0")

	if(dirname(output_pdb)!=""):
		submit_subprocess("rm %s " %(rosetta_ready_pdb) ) #neccesary if output_pdb is abspath..the rosetta_ready_pdb wouldn't be in the storage_folder
	#################################################################################

if(dirname(output_pdb)==""):
	main_folder_output_pdb=os.path.abspath("../%s" %(basename(output_pdb)))

	if(exists("%s" %(main_folder_output_pdb) )):
		print "Warning %s already exist...removing!" %(main_folder_output_pdb)
		submit_subprocess( "rm %s " %(main_folder_output_pdb) )	

	submit_subprocess("cp %s %s" %(output_pdb, main_folder_output_pdb) )

if(remove_storage_folder==True):
	os.chdir( main_folder )
	submit_subprocess( "rm -r %s " %(combine_pdb_folder) )	

print "SUCCESSFULLY COMPLETED!"

#pdb_segment={}

#parse the data

#print argv




'''
load /Users/sripakpa/minirosetta/Rosetta_rna_file/Square_RNA/idealize_helices/create_VDW_rep_screener/CG_BP.pdb
load /Users/sripakpa/minirosetta/Rosetta_rna_file/Square_RNA/idealize_helices/create_VDW_rep_screener/AU_BP.pdb

select base, name c2+c4+c5+c6+c8+n1+n2+n3+n4+n6+n7+n9+o2+o4+o6+n1p

align AU_BP and res 2 and base, short_corner_AB_1_10 and res 6 and base
align CG_BP and res 1 and base, short_corner_AB_1_10 and res 1 and base

SAVE ALIGN PDB (by mouse)

Fixed/Format the outputted aligned BP :  fix_pymol_output_pdb.py -input_pdb_file AU_BP_aligned_to_AB_corner.pdb

sripakpa@rescomp-09-171908:~/minirosetta/Formatted_RNA_pdb/Square_RNA/seperate_into_corners/short$ mv o2star_pack_short_corner_AB_1_10.pdb  short_corner_AB_1_10.pdb 
sripakpa@rescomp-09-171908:~/minirosetta/Formatted_RNA_pdb/Square_RNA/seperate_into_corners/short$ mv o2star_pack_short_corner_CD_11_20.pdb  short_corner_CD_11_20.pdb 
sripakpa@rescomp-09-171908:~/minirosetta/Formatted_RNA_pdb/Square_RNA/seperate_into_corners/short$ mv o2star_pack_short_corner_EF_21_30.pdb  short_corner_EF_21_30.pdb 
sripakpa@rescomp-09-171908:~/minirosetta/Formatted_RNA_pdb/Square_RNA/seperate_into_corners/short$ mv o2star_pack_short_corner_GH_31_40.pdb  short_corner_GH_31_40.pdb 
'''

'''
load /Users/sripakpa/minirosetta/Rosetta_rna_file/Square_RNA/create_VDW_rep_screener/C4G5_lower_VDW_rep_screener.pdb
load /Users/sripakpa/minirosetta/Rosetta_rna_file/Square_RNA/create_VDW_rep_screener/A4U5_upper_VDW_rep_screener.pdb

select base, name c2+c4+c5+c6+c8+n1+n2+n3+n4+n6+n7+n9+o2+o4+o6+n1p

align A4U5_upper_VDW_rep_screener and res 5 and base, START_AB_corner and res 8 and base
align C4G5_lower_VDW_rep_screener and res 4 and base, START_AB_corner and res 3 and base


'''
#combine_pdb.py AU_BP_aligned_to_AB_corner.pdb-1-1 CG_BP_aligned_to_CG_corner.pdb-2-2 short_corner_AB_1_10.pdb-1-6

