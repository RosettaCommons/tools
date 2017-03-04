#!/usr/bin/env python

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
import os
from os.path import basename, dirname, exists, expanduser
from time import sleep
import copy
import string
from os import popen 
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################

##############################################################


if(use_new_src_code()): 
	rna_swa_test_exe = get_rosetta_EXE("swa_rna_main") 
else:
	rna_swa_test_exe = get_rosetta_EXE("rna_swa_test") 

database_folder= get_rosetta_database_folder()

######################################################################



#calculate_RMSD.py -rmsd_res 16-23 -native_pdb ../3p59_RNA_A.pdb -alignment_res 16 -silent_file ../Nano_corner_SWA_models.out

#calculate_RMSD.py -rmsd_res 8-12 33-37 57-62 83-85 -native_pdb ../3p59_RNA_A.pdb -alignment_res 16 -silent_file ../Nano_corner_SWA_models.out


#calculate_RMSD.py -rmsd_res 1-15 26-40 51-65 76-88 -native_pdb ../3p59_RNA_A.pdb -alignment_res 16 -silent_file ../Nano_corner_SWA_models.out

#calculate_RMSD.py -rmsd_res 17-25 42-50 67-75 90-97 -native_pdb ../3p59_RNA_A.pdb -alignment_res 16 -silent_file ../Nano_corner_SWA_models.out

#calculate_RMSD.py -rmsd_res 1-97 -native_pdb ../3p59_RNA_A.pdb -alignment_res 16 -silent_file ../Nano_corner_SWA_models.out


#calculate_RMSD.py -rmsd_res 5-15 -native_pdb ../3p59_RNA_A.pdb -alignment_res 16 -silent_file ../Nano_corner_SWA_models.out


#calculate_RMSD.py -rmsd_res 8-12 -native_pdb ../3p59_RNA_A.pdb -alignment_res 16 -silent_file ../Nano_corner_SWA_models.out

#calculate_RMSD.py -rmsd_res 8 9 10 12 -native_pdb ../3p59_RNA_A.pdb -alignment_res 16 -silent_file ../Nano_corner_SWA_models.out

#calculate_RMSD.py -rmsd_res 12 13 -native_pdb ../3p59_RNA_A.pdb -alignment_res 16 -silent_file ../Nano_corner_SWA_models.out


#calculate_RMSD.py -rmsd_res 33-37 -native_pdb 3p59_RNA_A.pdb -alignment_res 16 -silent_file Nano_corner_SWA_models.out

#calculate_RMSD.py -rmsd_res 57-62 -native_pdb 3p59_RNA_A.pdb -alignment_res 16 -silent_file Nano_corner_SWA_models.out

#calculate_RMSD.py -rmsd_res 83-85 -native_pdb 3p59_RNA_A.pdb -alignment_res 16 -silent_file Nano_corner_SWA_models.out



rmsd_res_list = parse_segment_string_list(parse_options( argv, "rmsd_res", [""] ), verbose=True)

alignment_res_list =parse_options( argv, "alignment_res", [-1] )

if(len(rmsd_res_list)==0): error_exit_with_message("len(rmsd_res_list)==0")

if(len(alignment_res_list)==0): 
	alignment_res_list=[1]
	print "User did not specified align_res_list, using DEFUALT align_res_list=%s" %(alignment_res_list)


native_pdb = parse_options( argv, "native_pdb", "" )

if(native_pdb==""): error_exit_with_message('native_pdb==""')

if(exists(native_pdb)==False): error_exit_with_message("native_pdb (%s) doesn't exist!" %(native_pdb))


models_silent_file = parse_options( argv, "silent_file", "" )

if(exists(models_silent_file)==False): error_exit_with_message("models_silent_file (%s) doesn't exist!" %(models_silent_file) )


fasta_file="CALC_RMSD_fasta"
submit_subprocess('pdb2fasta.py %s > %s' %(native_pdb,fasta_file))

sequence = open( fasta_file ).readlines()[1][:-1] 

total_res=len(sequence)

all_res_list=range(1, total_res+1)

mock_sample_res=total_res

VERBOSE=True

output_filename="RECAL_RMSD_" + basename(models_silent_file)

if(exists(output_filename)): error_exit_with_message("output_filename already exist!")

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_args=%s " %(list_to_string(argv) ) )


command=rna_swa_test_exe
command += " -algorithm cluster_old"
command += " -in:file:silent " + models_silent_file
command += " -in:file:silent_struct_type  binary_rna"
command += " -database %s " %(database_folder)
command += " -out:file:silent " + output_filename
command += " -skip_clustering true "
command += " -clusterer_keep_pose_in_memory false " #Dec 11, 2011
command += " -clusterer_rename_tags false "
command += " -align_only_over_base_atoms true " 
command += " -simple_append_map true "
if(VERBOSE):
	command += " -VERBOSE true "
command += " -recreate_silent_struct true "

command += ' -native %s ' %(native_pdb)

command += " -rmsd_res %s " %(list_to_string(rmsd_res_list))
command += " -fixed_res %s " %(list_to_string( list(set(all_res_list)-set([mock_sample_res])) ) )
#command += " -fixed_res 1 "
command += " -input_res %s " %(list_to_string(all_res_list))
command += " -alignment_res %s " %(list_to_string(alignment_res_list, "-")[1:])
command += " -fasta %s " %(fasta_file)
command += " -sample_res %d " %(mock_sample_res)

command += " > CALC_RMSD_output.txt "

print "CALC_RMSD commnad=%s" %(command)

submit_subprocess(command)

submit_subprocess('grep "SCORE: " %s > SCORE_%s' %(output_filename, output_filename) )

score_lines = open( "SCORE_%s" %(output_filename), 'r').readlines()

description_line=score_lines[0]

score_lines=score_lines[1:]

score_lines=sorted(score_lines, key=lambda score_line: score_line.split()[10])

SCORE_LINES = open( "SCORE_%s" %(output_filename), 'w')

SCORE_LINES.write(description_line)

for line in score_lines:
	SCORE_LINES.write(line)

SCORE_LINES.close()	



