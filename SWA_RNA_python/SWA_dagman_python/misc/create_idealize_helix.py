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
######################################################################
print "---------FEB 03, 2012 WARNING: ...THIS FUNCTION IS DEPRECATED---------"
print "---------FEB 03, 2012 WARNING: ...THIS FUNCTION IS DEPRECATED---------"


#It seems that if you want to align a upper element of one motif with the lower element of the next motif, then old code lead to prefect alignment
#NEW CODE ACTUALLY MAKE ALIGNMENT WORST. SO for VC_II, use old code!

#The only situation that new_code make alignment better is if you have a symmetric duplex and we want prefect alignment between the lower and upper element..(this is basically inverting two times)
#A related point is that the current implementation of prefect_BP_symmetry can be used to exposed the assymmetry caused by the the mutation due to the invertion of the helix.

#OK, FINALLY realized that the input idealized helix (4_AU_BP.pdb) is not symmetric with respect to invertion.

#BOTH THESE CASE STILL DOESN'T SOLVE the short Hermann duplex problem. In that case you need a helix that is symmetric with respect to invertion.
#See ~/minirosetta/Rosetta_rna_file/VC_riboswitch/create_symmetric_helix_test/NO_mutation_Idealize_BP_alignment_test for proof.

#create_idealize_helix.py -fasta fasta -cutpoint_open 6 -prefect_BP_symmetry_mode True
#create_idealize_helix.py -fasta fasta -prefect_BP_symmetry_mode True

print "----------------------------------------------------------------------------------------------------------------------------------"
print "----------------------------------------------------------------------------------------------------------------------------------"
print "Enter: %s " %(list_to_string(argv))
print "----------------------------------------------------------------------------------------------------------------------------------"
print "----------------------------------------------------------------------------------------------------------------------------------"



HOMEDIR = expanduser('~') 

fasta_file = parse_options( argv, "fasta", "fasta" )
cutpoint_open= parse_options( argv, "cutpoint_open", 0 )
verbose= parse_options( argv, "verbose", "True" )
invert_helix= parse_options( argv, "invert_helix", "False" )

 
#if(prefect_BP_symmetry_mode==True): #Basically want a helix that is symmetric with itself.
#	error_exit_with_message("prefect_BP_symmetry_mode==True, this mode is not CORRECTLY implemented yet. Implementation below is wrong") 
	

four_AU_BP_pdb = get_PYPATH("database/IDEALIZE_HELIX/4_AU_BP.pdb") #Feb 03, 2012: This is deprecated 

four_BP_5_UG_3_wobble_pdb = get_PYPATH("database/IDEALIZE_HELIX/GENERIC_5_U_cutpoint_G_3_WOBBLE/4BP_5_UG_3_wobble.pdb")

invert_four_BP_5_UG_3_wobble_pdb = get_PYPATH("database/IDEALIZE_HELIX/GENERIC_5_U_cutpoint_G_3_WOBBLE/invert_chain_seq_cutpoint_4_4BP_5_UG_3_wobble.pdb")


###"WARNING: Temporarily using 5'-GU-3' wobble helical stub extract from 2R8S_J5_5A..until have time to generate generic version" 
four_BP_5_GU_3_wobble_pdb = get_PYPATH("database/IDEALIZE_HELIX/2R8S_J5_J5A_3_U_cutpoint_G_5_WOBBLE/2R8S_J5_JA_4BP_5_GU_3_wobble.pdb")
invert_four_BP_5_GU_3_wobble_pdb = get_PYPATH("database/IDEALIZE_HELIX/2R8S_J5_J5A_3_U_cutpoint_G_5_WOBBLE/invert_chain_seq_cutpoint_4_2R8S_J5_JA_4BP_5_GU_3_wobble.pdb")


storage_folder="mutate_helix/"

if( exists(storage_folder) ):
	print "Warning storage_folder (%s) already exist ... removing " %(storage_folder) 
	submit_subprocess("rm -r %s " %(storage_folder) )

submit_subprocess("mkdir %s " %(storage_folder) )

if(cutpoint_open==0):
	helix_type="lower"
else:
	helix_type="upper"

if(exists(fasta_file)==False): error_exit_with_message("fasta_file (%s) doesn't exist!!" %(fasta_file)) 

sequence = open( fasta_file  ).readlines()[1][:-1]
SEQUENCE=sequence.upper()

if(verbose): print "sequence= %s " %(sequence)
if(verbose): print "SEQUENCE= %s " %(SEQUENCE)

total_res=len(sequence)

if(verbose): print "total_res= %d " %(total_res) 

idealize_helix_pdb=""

res_1=0
res_2=0
res_3=0
res_4=0
mutate_res_string=""
element_slice_segment=""
VDW_screener_slice_segment=""
element_pdb_name="%s_elements.pdb" %(helix_type)
VDW_screener_pdb_name="%s_VDW_rep_screener.pdb" %(helix_type)

if(helix_type=="lower"):
	res_1=1
	res_2=2
	res_3=total_res-1
	res_4=total_res
	element_slice_segment="3 6"
	VDW_screener_slice_segment="1 3 6 8"

	if(invert_helix==False): #this is the original version of the code.

		if(SEQUENCE[res_2-1]=="U" and SEQUENCE[res_3-1]=="G"):
			idealize_helix_pdb=four_BP_5_UG_3_wobble_pdb
			mutate_res_string="%s3-%s6" %(SEQUENCE[res_1-1], SEQUENCE[res_4-1])
		elif(SEQUENCE[res_2-1]=="G" and SEQUENCE[res_3-1]=="U"):
			print "WARNING: Temporarily using 5'-GU-3' wobble helical stub extract from 2R8S_J5_5A..until have time to generate generic version" 
			idealize_helix_pdb=four_BP_5_GU_3_wobble_pdb
			mutate_res_string="%s3-%s6" %(SEQUENCE[res_1-1], SEQUENCE[res_4-1])
		else:
			idealize_helix_pdb=four_AU_BP_pdb
			mutate_res_string="%s3-%s4-%s5-%s6" %(SEQUENCE[res_1-1], SEQUENCE[res_2-1], SEQUENCE[res_3-1], SEQUENCE[res_4-1])

	else: #invert_helix

		if(SEQUENCE[res_2-1]=="U" or SEQUENCE[res_3-1]=="U"): #2nd BP is a 5'-GU-3' or 5'-UG-3' BP
			error_exit_with_message("invert_helix==True but GU wobble pairs..invert_helix is identical to standard_mode in this case")
		else:
			idealize_helix_pdb=four_AU_BP_pdb
			mutate_res_string="%s7-%s8-%s1-%s2" %(SEQUENCE[res_1-1], SEQUENCE[res_2-1], SEQUENCE[res_3-1], SEQUENCE[res_4-1])


else: #helix_type=="upper"

	res_1=cutpoint_open-1
	res_2=cutpoint_open
	res_3=cutpoint_open+1
	res_4=cutpoint_open+2

	element_slice_segment="1 2 7 8"
	VDW_screener_slice_segment="2 7"

	if(invert_helix==False): #this is the original version of the code.

		if(SEQUENCE[res_4-1]=="U" and SEQUENCE[res_1-1]=="G"):
			idealize_helix_pdb=invert_four_BP_5_UG_3_wobble_pdb
			mutate_res_string="%s2-%s7" %(SEQUENCE[res_2-1], SEQUENCE[res_3-1])
		elif(SEQUENCE[res_4-1]=="G" and SEQUENCE[res_1-1]=="U"):
			print "WARNING: Temporarily using 5'-GU-3' wobble helical stub extract from 2R8S_J5_5A..until have time to generate generic version" 
			idealize_helix_pdb=invert_four_BP_5_GU_3_wobble_pdb
			mutate_res_string="%s2-%s7" %(SEQUENCE[res_2-1], SEQUENCE[res_3-1])
		else:
			idealize_helix_pdb=four_AU_BP_pdb
			mutate_res_string="%s1-%s2-%s7-%s8" %(SEQUENCE[res_1-1], SEQUENCE[res_2-1], SEQUENCE[res_3-1], SEQUENCE[res_4-1])

	else: #invert_helix

		if(SEQUENCE[res_4-1]=="U" or SEQUENCE[res_1-1]=="U"):
			error_exit_with_message("invert_helix==True but GU wobble pairs..invert_helix is identical to standard_mode in this case")
		else:
			idealize_helix_pdb=four_AU_BP_pdb
			mutate_res_string="%s5-%s6-%s3-%s4" %(SEQUENCE[res_1-1], SEQUENCE[res_2-1], SEQUENCE[res_3-1], SEQUENCE[res_4-1])


		

#############################################################################################################################


if(idealize_helix_pdb==""): error_exit_with_message("idealize_helix_pdb==\"\"")

if(PATH_exists(idealize_helix_pdb)==False): error_exit_with_message("idealize_helix_pdb (%s) doesn't exist!!" %(idealize_helix_pdb))

if(idealize_helix_pdb==basename(idealize_helix_pdb)):
	error_exit_with_message("idealize_helix_pdb==basename(idealize_helix_pdb), idealize_helix_pdb=%s" %(idealize_helix_pdb))

submit_subprocess("cp %s %s " %(idealize_helix_pdb, basename(idealize_helix_pdb)) )

idealize_helix_pdb=basename(idealize_helix_pdb)

submit_subprocess("replace_chain_inplace.py  %s A " %(idealize_helix_pdb) )
#############################################################################################################################

SWA_mutate_command_line="SWA_mutate_residues.py " + idealize_helix_pdb + " " + mutate_res_string

if(verbose): print "SWA_mutate_command_line %s " %(SWA_mutate_command_line)

submit_subprocess(SWA_mutate_command_line)

mutate_helix_pdb="mutate_" + mutate_res_string + "_" + idealize_helix_pdb

if(exists(mutate_helix_pdb)==False):  error_exit_with_message("mutate_helix_pdb (%s) doesn't exist!!" %(mutate_helix_pdb))


#######################################################
if(invert_helix==True): #Nov 23, 2010. 
	
	inverted_mutate_helix_pdb="inverted_" + mutate_helix_pdb
	invert_helix_command="invert_chain_sequence.py -cutpoint_open %s -input_pdb_file %s -output_pdb_file %s " %(4, mutate_helix_pdb, inverted_mutate_helix_pdb) 
	#4 is the cutpoint open of the idealized helix and not the actual motif we are trying model (the one corresponding to the fasta file)

	if(verbose): print "invert_helix_command= %s" %(invert_helix_command)
	submit_subprocess(invert_helix_command)

	submit_subprocess("mv %s  %s/%s" %(mutate_helix_pdb,storage_folder, mutate_helix_pdb))
	mutate_helix_pdb=inverted_mutate_helix_pdb

#######################################################

create_element_command="SWA_pdbslice.py %s -segments %s %s" %(mutate_helix_pdb, element_slice_segment, element_pdb_name)
if(verbose): print "create_element_command= %s" %(create_element_command) 
submit_subprocess(create_element_command)

assert(exists("%s" %(element_pdb_name) ) )
#submit_subprocess("mv %s%s %s" %(element_pdb_name , mutate_helix_pdb, element_pdb_name) ) 

submit_subprocess("renumber_pdb_in_place.py %s" %( element_pdb_name) ) 

#######################################################

create_VDW_screener_command="SWA_pdbslice.py %s -segments %s %s" %(mutate_helix_pdb, VDW_screener_slice_segment, VDW_screener_pdb_name)
if(verbose): print "create_VDW_screener_command= %s" %(create_VDW_screener_command) 
submit_subprocess(create_VDW_screener_command)

assert(exists("%s" %(VDW_screener_pdb_name) ) )
#submit_subprocess("mv %s%s %s" %(VDW_screener_pdb_name, mutate_helix_pdb, VDW_screener_pdb_name) ) 

submit_subprocess("renumber_pdb_in_place.py %s" %( VDW_screener_pdb_name) ) 

#lower_VDW_screener.pdb  3-4 1-16 upper_VDW_screener.pdb  1-6 8-9 '
VDW_rep_screen_info_string= "-VDW_rep_screen_info %s" %(VDW_screener_pdb_name)
if(helix_type=="lower"):
	VDW_rep_screen_info_string+=" 3-4 %d-%d " %(res_1, res_4)
else:
	VDW_rep_screen_info_string+=" 1-6 %d-%d " %(res_2, res_3)

submit_subprocess(' echo "%s" >> VDW_rep_screen_info.txt' %(VDW_rep_screen_info_string) )
#######################################################


submit_subprocess("mv %s  %s/%s" %(mutate_helix_pdb,storage_folder, mutate_helix_pdb))
submit_subprocess("mv %s  %s/%s" %(idealize_helix_pdb,storage_folder, idealize_helix_pdb))

submit_subprocess("mv mutate_%s_%s_output.txt mutate_helix/mutate_%s_%s_output.txt" %(idealize_helix_pdb,mutate_res_string, idealize_helix_pdb,mutate_res_string))




