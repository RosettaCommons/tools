#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
import os
import copy


from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
from SWA_dagman_python.utility.RNA_BP_and_BS_util import *
from SWA_dagman_python.utility.extract_FR3D_data_functions import extract_FR3D_data_func

#extract_FR3D_data.py -pdb_name  1y26_RNA_A.pdb  -sample_segment 1-71 -double_count_base_pair false -INCLUDE_EDGE_PHOSPHATE false


#extract_FR3D_data.py -pdb_name  S_000004.pdb   -sample_segment 5-9 > LOG_extract_FR3D_data.txt
#extract_FR3D_data.py -pdb_name  S_000000.pdb -sample_segment 22-28> LOG_extract_FR3D_data.txt
#extract_FR3D_data.py -pdb_name  S_000004.pdb  -sample_segment 3-9 > LOG_extract_FR3D_data.txt
#extract_FR3D_data.py -pdb_name  S_000001.pdb -sample_segment 50-59 > LOG_extract_FR3D_data.txt


#extract_FR3D_data.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/LOOP_2PN4_7_11/expand_radius_50_2pn4_RNA_A.pdb   -sample_segment 3-7 > LOG_extract_FR3D_data.txt

#extract_FR3D_data.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/LOOP_1NUJ_54_57/expand_radius_50_1nuj_RNA_A.pdb   -sample_segment 8-11 > LOG_extract_FR3D_data.txt
#extract_FR3D_data.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/LOOP_2QBZ_142_145/expand_radius_50_2qbz_RNA_A.pdb  -sample_segment 10-13 > LOG_extract_FR3D_data.txt

#extract_FR3D_data.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/LOOP_1U6B_203_207/expand_radius_50_1u6b_RNA_A.pdb   -sample_segment 8-12 > LOG_extract_FR3D_data.txt
#extract_FR3D_data.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/LOOP_2QWY_60_66/expand_radius_50_combine_2qwy.pdb  -sample_segment 3-9 > LOG_extract_FR3D_data.txt
#extract_FR3D_data.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/LOOP_3D2V_117_121/expand_radius_50_3d2v_RNA_A.pdb  -sample_segment 10-14 > LOG_extract_FR3D_data.txt

#extract_FR3D_data.py -pdb_name  Feb_12_SWA_LOOP_3D2V_117_121_S_0_000001.pdb  -sample_segment 10-14 > LOG_extract_FR3D_data.txt


#extract_FR3D_data.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/LOOP_3OWI_161_167/expand_radius_50_3owi_RNA_A.pdb -sample_segment 22-28 > LOG_extract_FR3D_data.txt
#extract_FR3D_data.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/LOOP_3G78_286_292/expand_radius_50_3g78_RNA_A.pdb  -sample_segment 4-10 > LOG_extract_FR3D_data.txt

#extract_FR3D_data.py -pdb_name  LOOP_3G78_286_292_S_000002.pdb  -sample_segment 4-10 > LOG_extract_FR3D_data.txt

#extract_FR3D_data.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/S1_2R8S_J5_5a/J5_5a_2r8s.pdb   -sample_segment 5-9 > LOG_extract_FR3D_data.txt

#extract_FR3D_data.py -pdb_name  L1_2R8S_S_4.pdb  -sample_segment 5-9 > LOG_extract_FR3D_data.txt

#extract_FR3D_data.py -pdb_name  J5_J5a_working_native_pose.pdb   -sample_segment 5-9 > LOG_extract_FR3D_data.txt

#extract_FR3D_data.py -pdb_name  Mar_J5_Ja_S_3_5.pdb -sample_segment 5-9 > LOG_extract_FR3D_data.txt



 
#extract_FR3D_data.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/S2_2R8S_J5_5a/J5_5a_2r8s.pdb   -sample_segment 14-17 > LOG_extract_FR3D_data.txt
#extract_FR3D_data.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/1L2X/1l2x_RNA_B_with_sym_partner.pdb  -sample_segment 18-24 > LOG_extract_FR3D_data.txt
#extract_FR3D_data.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/1S72_23_rRNA_35_40/expand_radius_100_1s72_RNA_A_35_40_3_8.pdb   -sample_segment 3-8 > LOG_extract_FR3D_data.txt
#extract_FR3D_data.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/1S72_23_rRNA_520_525/expand_radius_100_1s72_RNA_A_520_525_10_15.pdb   -sample_segment 10-15 > LOG_extract_FR3D_data.txt
#extract_FR3D_data.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/1S72_23_rRNA_1950_1959/expand_radius_100_1s72_RNA_A_1950_1959_50_59.pdb   -sample_segment 50-59 > LOG_extract_FR3D_data.txt
#extract_FR3D_data.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/1S72_23_rRNA_1923_1932/expand_radius_100_1s72_RNA_A_1923_1932_15_24.pdb   -sample_segment 15-24 > LOG_extract_FR3D_data.txt
#extract_FR3D_data.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/1S72_23_rRNA_2376_2382/expand_radius_100_1s72_RNA_A_2376_2382_72_78.pdb   -sample_segment 72-78 > LOG_extract_FR3D_data.txt

#extract_FR3D_data.py -pdb_name 1s72_RNA_A_1950_1959_50_59_R_000001.pdb   -sample_segment 50-59

#double_count_base_pair

 


#extract_FR3D_data.py -pdb_name expand_radius_100_1s72_RNA_A_520_525_10_15.pdb -sample_segment 10-15
#extract_FR3D_data.py -pdb_name expand_radius_100_1s72_RNA_A_1950_1959_50_59.pdb  -sample_segment 50-59

#extract_FR3D_data.py -pdb_name 10bp_polyua_aform_rna_generated_by_3dnaprogram_RNA_A.pdb -sample_segment 2-19
#extract_FR3D_data.py -pdb_name mutate_C1-C2-C3-C4-C5-C6-C7-C8-C9-C10-G11-G12-G13-G14-G15-G16-G17-G18-G19-G20_10bp_polyua_aform_rna_generated_by_3dnaprogram_RNA_A.pdb -sample_segment 2-19

#extract_FR3D_data.py  -pdb_name mutate_U15-A16_SWA_best_energy.pdb  -sample_segment All -double_count False

#extract_FR3D_data.py -pdb_name mutate_U15-A16_FARFAR_S_0_38.pdb  -sample_segment All -double_count False



copy_argv=copy.deepcopy(argv)

HOMEDIR = expanduser('~') 

print "HOMEDIR= %s" %(HOMEDIR)


#####################################

sample_segment= parse_options( argv, "sample_segment", "" )

run_rosetta_extract_hbond= parse_options( argv, "run_rosetta_extract_hbond", "False" )

rosetta_pdb_name= parse_options( argv, "pdb_name", "1s72_RNA_A.pdb" )

double_count= parse_options( argv, "double_count", "True" ) #Use True for single stranded loops (Aug 11, 2011)

if(run_rosetta_extract_hbond):
	if(double_count==False): error_exit_with_message("double_count=False mode not yet implemented in run_rosetta_extract_hbond!")
	submit_subprocess("extract_hbond_stat.py -pdb_name %s -sample_segment %s " %(rosetta_pdb_name, sample_segment ) )

if( dirname(rosetta_pdb_name)!="" ):rosetta_pdb_name=HOMEDIR + rosetta_pdb_name

allow_near_match= parse_options( argv, "allow_near_match", "False" )

extract_FR3D_data_func(rosetta_pdb_name, sample_segment, double_count, allow_near_match, Verbose=True)

print 

print "----------------------------------------------------------------------------------------------------------------------------"
print "----------------------%s----------------------" %(list_to_string(copy_argv))
print "extract_FR3D_data.py sucessfully RAN! "
print "----------------------------------------------------------------------------------------------------------------------------"



