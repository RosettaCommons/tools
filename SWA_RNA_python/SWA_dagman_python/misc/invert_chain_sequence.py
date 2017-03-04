#!/usr/bin/env python

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
from os.path import basename, dirname, exists, expanduser
from time import sleep


from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options

#invert_chain_sequence.py -cutpoint_open 5 -input_pdb_file  1f5h_upper_duplex.pdb

#invert_chain_sequence.py -cutpoint_open 8 -input_pdb_file  A_site_1t0e_RNA_B.pdb

#invert_chain_sequence.py -cutpoint_open 23 -input_pdb_file  daslab_ts001_1_RNA_A.pdb 

#invert_chain_sequence.py -cutpoint_open 4 -input_pdb_file  2R8S_J5_JA_4BP_5_GU_3_wobble.pdb 


#invert_chain_sequence.py -cutpoint_open 11 -input_pdb_file 2PN3_Square_RNA.pdb

#invert_chain_sequence.py -cutpoint_open 10 -input_pdb_file 10bp_polyua_aform_rna_generated_by_3dnaprogram_RNA_A.pdb
#invert_chain_sequence.py -cutpoint_open 2 -input_pdb_file  upper_element.pdb
#invert_chain_sequence.py -cutpoint_open 7 -input_pdb_file  upper_VDW_rep_screen.pdb 

HOMEDIR = expanduser('~') 

cutpoint_open= parse_options( argv, "cutpoint_open", 0 )
input_pdb_file= parse_options( argv, "input_pdb_file", "" )

output_pdb_file= parse_options( argv, "output_pdb_file", "")

if(cutpoint_open==0): error_exit_with_message("cutpoint_open==0") 

if(input_pdb_file==""): error_exit_with_message("input_pdb_file==\"\"") 


if(exists(input_pdb_file)==False): error_exit_with_message("input_pdb_file (%s) doesn't exist!!" %(input_pdb_file))


local_input_pdb_file="local_invert_chain_sequence_input_pdb_file.pdb"

if(exists(local_input_pdb_file)==True): 
	print "WARNING: local_input_pdb_file (%s) already exist!! ...removing..." %(local_input_pdb_file)
	submit_subprocess("rm %s " %(local_input_pdb_file) )

submit_subprocess("cp %s %s " %(input_pdb_file,local_input_pdb_file) )

if(output_pdb_file==""):
	output_pdb_file="invert_chain_seq_cutpoint_%s_%s" %(cutpoint_open, basename(input_pdb_file))


input_pdb_file=""

submit_subprocess("renumber_pdb_in_place.py %s" %( local_input_pdb_file) ) 

fasta_file="fasta_" + local_input_pdb_file.replace(".pdb","")

submit_subprocess('SWA_pdb2fasta.py %s > %s' %(local_input_pdb_file,fasta_file))

sequence = open( fasta_file ).readlines()[1][:-1]
#if(verbose): print "sequence= %s " %(sequence)

total_res=len(sequence)
print "total_res= %s" %(total_res)


new_to_old_seq_num_map={}

new_cutpoint_open=total_res-cutpoint_open

print "orig_cutpoint_open= %s " %(cutpoint_open)
print "new_cutpoint_open= %s " %(new_cutpoint_open)


for new_seq_num in range(1,total_res+1):

	if(new_seq_num<=new_cutpoint_open):
		old_seq_num=new_seq_num+cutpoint_open
	else: 
		old_seq_num=new_seq_num-new_cutpoint_open

	print "new_seq_num=%s, old_seq_num=%s" %(new_seq_num, old_seq_num)
	new_to_old_seq_num_map[new_seq_num]=old_seq_num


print "new_to_old_seq_num_map: ", new_to_old_seq_num_map


if(exists(output_pdb_file)==True): 
	print "WARNING: output_pdb_file (%s) already exist!! ...removing..." %(output_pdb_file)
	submit_subprocess("rm %s " %(output_pdb_file) )

OUTPUT_PDB = open(output_pdb_file, 'w')

input_pdb_line_list = open( local_input_pdb_file ).readlines()

for new_seq_num in range(1,total_res+1):

	found_pdb_line=False

	for input_pdb_line in input_pdb_line_list:
	
		input_pdb_line_split=input_pdb_line.split()
		pdb_line_seq_num=int(input_pdb_line_split[5])

		if(new_seq_num not in new_to_old_seq_num_map): error_exit_with_message("new_seq_num (%s) is not in new_to_old_seq_num_map!!" %(new_seq_num))
		old_seq_num=new_to_old_seq_num_map[new_seq_num]

		if(pdb_line_seq_num==old_seq_num):
			found_pdb_line=True
			OUTPUT_PDB.write( "%s\n" %(input_pdb_line ) )

		#print input_pdb_line_split, " line_seq_num= ", line_seq_num
	if(found_pdb_line==False): error_exit_with_message("found_pdb_line==False for new_seq_num= %s, old_seq_num= %s !!" %(new_seq_num, old_seq_num ) )
OUTPUT_PDB.close()

submit_subprocess("renumber_pdb_in_place.py %s" %(output_pdb_file) ) 


