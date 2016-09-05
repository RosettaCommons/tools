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


#extract_hbond_stat.py -pdb_name  1y26_RNA_A.pdb  -sample_segment 1-71 -double_count_base_pair false -INCLUDE_EDGE_PHOSPHATE false


#extract_hbond_stat.py -pdb_name  S_000004.pdb   -sample_segment 5-9
#extract_hbond_stat.py -pdb_name  S_000000.pdb -sample_segment 22-28
#extract_hbond_stat.py -pdb_name  S_000004.pdb  -sample_segment 3-9
#extract_hbond_stat.py -pdb_name  S_000001.pdb -sample_segment 50-59


#extract_hbond_stat.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/LOOP_2PN4_7_11/expand_radius_50_2pn4_RNA_A.pdb   -sample_segment 3-7

#extract_hbond_stat.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/LOOP_1NUJ_54_57/expand_radius_50_1nuj_RNA_A.pdb   -sample_segment 8-11
#extract_hbond_stat.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/LOOP_2QBZ_142_145/expand_radius_50_2qbz_RNA_A.pdb  -sample_segment 10-13

#extract_hbond_stat.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/LOOP_1U6B_203_207/expand_radius_50_1u6b_RNA_A.pdb   -sample_segment 8-12
#extract_hbond_stat.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/LOOP_2QWY_60_66/expand_radius_50_combine_2qwy.pdb  -sample_segment 3-9
#extract_hbond_stat.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/LOOP_3D2V_117_121/expand_radius_50_3d2v_RNA_A.pdb  -sample_segment 10-14
#extract_hbond_stat.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/LOOP_3OWI_161_167/expand_radius_50_3owi_RNA_A.pdb -sample_segment 22-28
#extract_hbond_stat.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/LOOP_3G78_286_292/expand_radius_50_3g78_RNA_A.pdb  -sample_segment 4-10
#extract_hbond_stat.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/S1_2R8S_J5_5a/J5_5a_2r8s.pdb   -sample_segment 5-9
#extract_hbond_stat.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/S2_2R8S_J5_5a/J5_5a_2r8s.pdb   -sample_segment 14-17
#extract_hbond_stat.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/1L2X/1l2x_RNA_B_with_sym_partner.pdb  -sample_segment 18-24
#extract_hbond_stat.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/1S72_23_rRNA_35_40/expand_radius_100_1s72_RNA_A_35_40_3_8.pdb   -sample_segment 3-8
#extract_hbond_stat.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/1S72_23_rRNA_520_525/expand_radius_100_1s72_RNA_A_520_525_10_15.pdb   -sample_segment 10-15
#extract_hbond_stat.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/1S72_23_rRNA_1950_1959/expand_radius_100_1s72_RNA_A_1950_1959_50_59.pdb   -sample_segment 50-59
#extract_hbond_stat.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/1S72_23_rRNA_1923_1932/expand_radius_100_1s72_RNA_A_1923_1932_15_24.pdb   -sample_segment 15-24
#extract_hbond_stat.py -pdb_name  /minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/1S72_23_rRNA_2376_2382/expand_radius_100_1s72_RNA_A_2376_2382_72_78.pdb   -sample_segment 72-78

#double_count_base_pair

#extract_hbond_stat.py -pdb_name expand_radius_100_1s72_RNA_A_520_525_10_15.pdb -sample_segment 10-15
#extract_hbond_stat.py -pdb_name expand_radius_100_1s72_RNA_A_1950_1959_50_59.pdb  -sample_segment 50-59

#extract_hbond_stat.py -pdb_name 10bp_polyua_aform_rna_generated_by_3dnaprogram_RNA_A.pdb -sample_segment 2-19
#extract_hbond_stat.py -pdb_name mutate_C1-C2-C3-C4-C5-C6-C7-C8-C9-C10-G11-G12-G13-G14-G15-G16-G17-G18-G19-G20_10bp_polyua_aform_rna_generated_by_3dnaprogram_RNA_A.pdb -sample_segment 2-19

#####################################################

copy_argv=copy.deepcopy(argv)

HOMEDIR = expanduser('~') 

print "HOMEDIR= %s" %(HOMEDIR)
#####################################################

EXE = get_rosetta_EXE("parin_test") 

database_folder= get_rosetta_database_folder()

####################################################

INCLUDE_EDGE_PHOSPHATE=parse_options( argv, "INCLUDE_EDGE_PHOSPHATE", "true" )

MODE=parse_options( argv, "MODE", "local") 

if(MODE!="manual_sub" and MODE!="standard_queue" and MODE!="local"):
	error_exit_with_message("Invalid MODE(%s)!" %(MODE))

sample_segment= parse_options( argv, "sample_segment", "" )

double_count_base_pair= parse_options( argv, "double_count_base_pair", "true" )

pdb_name= parse_options( argv, "pdb_name", "1s72_RNA_A.pdb" )

if(dirname(pdb_name)!=""): pdb_name=HOMEDIR + pdb_name

pdb_name=os.path.abspath(pdb_name)


if( exists( pdb_name )==False):  error_exit_with_message("pdb_name (%s) doesn't exist!" %(pdb_name) )

###########################consistency check#######################################

#start_res_list=[]
#end_res_list=[]


sample_segment=sample_segment.split('-')
if(len(sample_segment)!=2): error_exit_with_message("len(sample_segment)!=2" )

start_res=int(sample_segment[0])
end_res=int(sample_segment[1])

if(start_res<0):  error_exit_with_message("start_res is not a positive integer, start_res=%d" %(start_res) )
if(end_res<0):  error_exit_with_message("end_res is not positive integer, end_res=%d" %(end_res) )

if(start_res>=end_res): error_exit_with_message("start_res(%d) > end_res(%d)" %(start_res, end_res))

#start_res_list.append(start_res)
#end_res_list.append(end_res)
#if(len(start_res_list)!=len(end_res_list)):
#	print "start_res_list=", start_res_list
#	print "end_res_list=", end_res_list
#	error_exit_with_message("len(start_res_list)!=len(end_res_list)" )
##############################################################################


foldername="extract_hbond_%s_%4s_%4s" %( basename(pdb_name).replace(".pdb",""), str(start_res).zfill(4), str(end_res).zfill(4))

if( exists(foldername) ):   
	print "Warning foldername (%s) already exist! ...removing!..." %(foldername)
	submit_subprocess("rm -r %s " %(foldername) )

submit_subprocess( 'mkdir %s' %(foldername)  )

os.chdir( foldername)


sample_res=''

for seq_num in range(start_res, end_res+1):
	sample_res+='%d ' %(seq_num)


command=EXE
command += " -algorithm extract_hydrogen_bonds_statistic"
command += " -database %s" %(database_folder)
command += " -s " + pdb_name
command += " -sample_res " + sample_res
command += " -INCLUDE_EDGE_PHOSPHATE %s " %(INCLUDE_EDGE_PHOSPHATE) 
command += " -double_count_base_pair %s " %(double_count_base_pair)

job_tag = abspath('.').replace('/','_')

	
if(MODE=="manual_sub"):

	command ="bsub -q IA -Is  -o log.out -e log.err -J %s %s " %(job_tag,command)
	print command
	COMMAND= open( "COMMAND.txt", 'w')
	COMMAND.write(command + "\n")
	COMMAND.close()

elif(MODE=="standard_queue"):

	command ="bsub -W 140:0 -M 4000000 -o log.out -e log.err -J %s %s " %(job_tag,command)
	print command
	submit_subprocess( command )

elif(MODE=="local"):	

	command += " > log.out 2> log.err"
	print command
	submit_subprocess( command )

else:
	error_exit_with_message("Invalid mode=%s "%(MODE))

	
os.chdir( "../")

print 

print "----------------------------------------------------------------------------------------------------------------------------"
print "----------------------%s----------------------" %(list_to_string(copy_argv))
print "extract_hbond_stat.py sucessfully RAN! "
print "----------------------------------------------------------------------------------------------------------------------------"

