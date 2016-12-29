#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
import os
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################

native_src =  parse_options( argv, "native_src", "" )  #native_src is silent_file.out-tag_name
decoy_src =	  parse_options( argv, "decoy_src", [""] )

alignment_res_pairs =	  parse_options( argv, "alignment_res_pairs", [""] ) #1-5 7-9  (res 1 of native to res 5 of decoy, res 7 of native to res 9 of decoy)
rmsd_res_pairs =	  parse_options( argv, "rmsd_res_pairs", [""] )

align_only_over_base_atoms=parse_options( argv, "align_only_over_base_atoms", "true")
dump_pdb=parse_options(argv, "dump_pdb", "true")
alignment_RMSD_CUTOFF =	  parse_options( argv, "alignment_RMSD_CUTOFF", 0.1 ) #ensure  perfect alignment

 
###-ignore_virtual_res ??

#calculate_pairwise_RMSD.py -dump_pdb true

##############################################################

native_src="Sigel_39mer_NMR_model.out-Sigel_39mer_NMR_model"
decoy_src=["3X2_loop_2aht.out-3X2_loop_2aht"]#"M_1.out-M_1", 
rmsd_res_pairs=["2-2","3-3","4-4","5-5","8-8","9-9","10-10","11-11","12-12"] 
alignment_res_pairs=["2-2","3-3","4-4","5-5","8-8","9-9","10-10","11-11","12-12"]
alignment_RMSD_CUTOFF=99.99



##############################################################
'''
native_src="/Volumes/Jan_2011_Ext_HD_Parin/minirosetta/02_2011/Feb_13_SWA_LOOP_3G78_286_292/denovo/main_folder/LOOP_3G78_286_292/TOP_ENERGY_CLUSTERS/top_energy_clusters.out-S_000000"
decoy_src=["/Volumes/Jan_2011_Ext_HD_Parin/minirosetta/Aug_1_ELLIPSE_3G78_LOOP/main_folder/Ellipse_LOOP_3G78_286_292/TOP_ENERGY_CLUSTERS/top_energy_clusters.out-S_000000"]
rmsd_res_pairs=["4-64", "5-65", "6-66", "7-67", "8-68", "9-69", "10-70"] 
alignment_res_pairs=["1-4", "33-138"]
alignment_RMSD_CUTOFF=0.1

#alignment_res_pairs=["3-63", "4-64", "5-65", "6-66", "7-67", "8-68", "9-69", "10-70", "11-71"] 
#alignment_RMSD_CUTOFF=1.0

#3G78 ENDS ALIGN: RMSD between native_tag(S_000000) and decoy_tag(S_000000) IS 0.125051
#3G78 LOOP_PLUS_EDGE :RMSD between native_tag(S_000000) and decoy_tag(S_000000) IS 0.124333
'''
##############################################################
'''
native_src="/Volumes/Jan_2011_Ext_HD_Parin/minirosetta/02_2011/Feb_12_SWA_LOOP_3D2V_117_121/denovo/main_folder/LOOP_3D2V_117_121/TOP_ENERGY_CLUSTERS/top_energy_clusters.out-S_000003"
decoy_src=["/Volumes/Jan_2011_Ext_HD_Parin/minirosetta/July_31_ELLIPSE_3D2V_LOOP/main_folder/Ellipse_LOOP_3D2V_117_121/TOP_ENERGY_CLUSTERS/top_energy_clusters.out-S_000003"]
rmsd_res_pairs=["10-38", "11-39", "12-40", "13-41", "14-42"]
alignment_res_pairs=["1-10", "26-69"]
alignment_RMSD_CUTOFF=0.1

#alignment_res_pairs=["9-37", "10-38", "11-39", "12-40", "13-41", "14-42", "15-43"]
#3D2V ENDS ALIGN: RMSD between native_tag(S_000003) and decoy_tag(S_000003) IS 0.138334
#3D2V LOOP_PLUS_EDGE ALIGN: RMSD between native_tag(S_000003) and decoy_tag(S_000003) IS 0.116335 
'''

##############################################################
'''
native_src="/Volumes/Jan_2011_Ext_HD_Parin/minirosetta/03_2011/Mar_18_SWA_J5_5a_LOOP/main_folder/S1_2R8S_J5_5a/TOP_ENERGY_CLUSTERS/top_energy_clusters.out-S_000003"
decoy_src=["/Volumes/Jan_2011_Ext_HD_Parin/minirosetta/Aug_16_Ellipse_S1_2R8S_21_25/main_folder/Ellipse_S1_2R8S_21_25/TOP_ENERGY_CLUSTERS/top_energy_clusters.out-S_000005"]
rmsd_res_pairs=["5-14", "6-15", "7-16", "8-17", "9-18"]
alignment_res_pairs=["1-10","21-49"]
alignment_RMSD_CUTOFF=0.1

#alignment_res_pairs=["4-13","5-14", "6-15", "7-16", "8-17", "9-18", "10-19"]
#S1_2R8S_21_25 ENDS ALIGN: RMSD between native_tag(S_000003) and decoy_tag(S_000005) IS 0.227545
#S1_2R8S_21_25 LOOP_PLUS_EDGE ALIGN: RMSD between native_tag(S_000003) and decoy_tag(S_000005) IS 0.196078
'''
##############################################################
'''
native_src="/Volumes/Jan_2011_Ext_HD_Parin/minirosetta/02_2011/Feb_10_2011_SWA_C7_2_THREE_nucleotide_FIX_LOWER_BULGE/pose_cluster_200_region_FINAL.out/Virtualize_U19_S_000000.out-Virtual_U19_S_0"
decoy_src=["/Volumes/Jan_2011_Ext_HD_Parin/minirosetta/Aug_30_SWA_extended_C7_2/region_FINAL.out-S_0"]

#decoy_src=["/Volumes/Parin_Ext_Sept_2010/minirosetta/03_and_04_2010_Success_Tetraloop_Receptor_RUNS/Apr6_SUCCESS_C7_2_full/region_FINAL.out-S_0"]
rmsd_res_pairs=["10-10", "11-11", "12-12", "13-13", "18-18", "19-19"] 
alignment_res_pairs=["9-9","10-10", "11-11", "12-12", "13-13","14-14", "17-17", "18-18", "19-19","20-20"]
alignment_RMSD_CUTOFF=2.0


#rmsd_res_pairs=["11-11", "12-12", "13-13"] 
#alignment_res_pairs=["10-10", "11-11", "12-12", "13-13", "14-14"]
#alignment_RMSD_CUTOFF=1.0

#alignment_res_pairs=["1-1", "22-22"]
#alignment_RMSD_CUTOFF=0.1

#C72 LOOP_PLUS_EDGE ALIGN :RMSD between native_tag(S_0) and decoy_tag(S_0) IS 0.915682
#C72 ENDS ALIGN           :RMSD between native_tag(S_0) and decoy_tag(S_0) IS 1.95566
'''

'''
native_src="/Volumes/Jan_2011_Ext_HD_Parin/minirosetta/03_2011/Mar_28_DOWNLOAD_Mar_27_FARFAR_C7_2_THREE_nucleotide_FIX_LOWER_BULGE/TRAIL_2_BIOX_SUBMIT/pose_cluster_200_FINAL_filtered_energy.out/Virtualize_U19_S_000000.out-Virtual_U19_S_0"
decoy_src=["/Volumes/Jan_2011_Ext_HD_Parin/minirosetta/Aug_29_EXTENDED_C7_2_FARFAR/pose_cluster_200_FINAL_filtered_energy.out/cluster_200_FINAL_filtered_energy.out-S_000000"]

rmsd_res_pairs=["10-10", "11-11", "12-12", "13-13", "18-18", "19-19"] 
alignment_res_pairs=["9-9","10-10", "11-11", "12-12", "13-13","14-14", "17-17", "18-18", "19-19","20-20"]
alignment_RMSD_CUTOFF=4.0

#Energy gap between first and second FARFAR extended C7.2 cluster center:  (-124.86200) + 123.94700 = -0.915

#alignment_res_pairs=["10-10", "11-11", "12-12", "13-13", "14-14"]
#alignment_RMSD_CUTOFF=4.0
#decoy_src=["/Volumes/Jan_2011_Ext_HD_Parin/minirosetta/Aug_24_TRAIL_2_EXTENDED_C7_2_FARFAR/CC_GG_helix/pose_cluster_200_FINAL_filtered_energy.out/cluster_200_FINAL_filtered_energy.out-S_000000"]
'''

#######################P4-P6 J5_J5a junction compare 2R8S (chain R) and 1GID (chain A) [note that chain B is missing nts and hence ignored############################
'''
native_src=HOMEDIR + "/minirosetta/test/Aug_25_annotate_full_legth_P4_P6_and_TLR/FULL_LENGTH/Virtualize_A24_2R8S_R.out-Virtual_A24_2R8S_R"
decoy_src=[HOMEDIR + "/minirosetta/test/Aug_25_annotate_full_legth_P4_P6_and_TLR/FULL_LENGTH/Virtualize_A24_1GID_A.out-Virtual_A22_1GID_A"]
#alignment_res_pairs=["20-18", "21-19", "22-20", "23-21", "24-22", "25-23", "26-24"]
alignment_res_pairs=["27-25"]
rmsd_res_pairs=["21-19", "22-20", "23-21", "24-22", "25-23"]
alignment_RMSD_CUTOFF=2.0
'''


##############################################################
'''
native_src="../No_bulge_native.out-native"
decoy_src=["../No_bulge_das_1.out-das_1"]
alignment_res_pairs=["2-2","3-3", "4-4", "5-5", "6-6", "9-9", "10-10", "11-11", "12-12", "13-13", "14-14"]
rmsd_res_pairs=["2-2","3-3", "4-4", "5-5", "10-10", "11-11", "12-12", "13-13", "14-14"]
alignment_RMSD_CUTOFF=10.0
align_only_over_base_atoms=False
'''
##############################################################


if(alignment_res_pairs==[""]): error_exit_with_message("alignment_res_pairs==[\"\"]")
if(rmsd_res_pairs==[""]): error_exit_with_message("rmsd_res_pairs==[\"\"]")

####Consistency check############################
for n in range(len(alignment_res_pairs)):
	alignment_res_pair=alignment_res_pairs[n].split("-")
	if(len(alignment_res_pair)!=2): error_exit_with_message("len(alignment_res_pair)!=2!, alignment_res_pairs[%d]=%s" %(n, alignment_res_pairs[n]) )
	native_res=int(alignment_res_pair[0])
	decoy_res=int(alignment_res_pair[1])

for n in range(len(rmsd_res_pairs)):
	rmsd_res_pair=rmsd_res_pairs[n].split("-")
	if(len(rmsd_res_pair)!=2): error_exit_with_message("len(rmsd_res_pair)!=2!, rmsd_res_pairs[%d]=%s" %(n, rmsd_res_pairs[n]) )
	native_res=int(rmsd_res_pair[0])
	decoy_res=int(rmsd_res_pair[1])

##################################################



if(native_src==""): error_exit_with_message("User need to pass in native_src!")
if(decoy_src==[""]): error_exit_with_message("User need to pass in decoy_src")

if( len(native_src.split("-"))!=2 ): error_exit_with_message("len(native_src.split("-"))!=2") 
if( len(decoy_src[0].split("-"))!=2 ):  error_exit_with_message("len(decoy_src[0].split("-"))!=2") 

native_silent_file=native_src.split("-")[0]
native_tag=        native_src.split("-")[1]

decoy_silent_file= decoy_src[0].split("-")[0]
decoy_tag=         decoy_src[0].split("-")[1]

if(dirname(native_silent_file)!=""): native_silent_file=os.path.abspath(native_silent_file)
if(dirname(decoy_silent_file)!=""): decoy_silent_file=os.path.abspath(decoy_silent_file)


if(exists(native_silent_file)==False): error_exit_with_message("native_silent_file (%s) doesn't exist!" %(native_silent_file) )
if(exists(decoy_silent_file)==False): error_exit_with_message("decoy_silent_file (%s) doesn't exist!" %(decoy_silent_file) )

	
if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_args=%s " %(list_to_string(argv) ) )

##############################################################



if(use_new_src_code()): 
	EXE = get_rosetta_EXE("swa_rna_util") 
else:
	EXE = get_rosetta_EXE("parin_test") 

##############################################################

command=EXE
command += ' -algorithm calculate_pairwise_RMSD'
command += ' -database %s' %(get_rosetta_database_folder())
command += ' -native %s' %(native_silent_file)
command += ' -native_tag_name %s' %(native_tag)
command += ' -s %s' %(decoy_silent_file)
command += ' -decoy_tag_name %s' %(decoy_tag)
command += ' -alignment_res_pairs %s' %(list_to_string(alignment_res_pairs))
command += ' -alignment_RMSD_CUTOFF %s' %(alignment_RMSD_CUTOFF)
command += ' -rmsd_res_pairs %s' %(list_to_string(rmsd_res_pairs))

command += ' -align_only_over_base_atoms %s ' %(align_only_over_base_atoms)
command += ' -dump %s ' %(dump_pdb)
command += ' -output_virtual  true '
command += ' > calculate_pairwise_RMSD.txt '

print command 
submit_subprocess( command )




