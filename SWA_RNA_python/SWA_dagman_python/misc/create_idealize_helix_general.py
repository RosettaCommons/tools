#!/usr/bin/env python

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
import os
from os.path import basename, dirname, exists, expanduser
from time import sleep

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
from SWA_dagman_python.utility.RNA_sequence import Is_perfect_WC_helix
######################################################################


#It seems that if you want to align a upper element of one motif with the lower element of the next motif, then old code lead to prefect alignment
#NEW CODE ACTUALLY MAKE ALIGNMENT WORST. SO for VC_II, use old code!

#The only situation that new_code make alignment better is if you have a symmetric duplex and we want prefect alignment between the lower and upper element..(this is basically inverting two times)
#A related point is that the current implementation of prefect_BP_symmetry can be used to exposed the assymmetry caused by the the mutation due to the invertion of the helix.

#OK, FINALLY realized that the input idealized helix (4_AU_BP.pdb) is not symmetric with respect to invertion.

#BOTH THESE CASE STILL DOESN'T SOLVE the short Hermann duplex problem. In that case you need a helix that is symmetric with respect to invertion.
#See ~/minirosetta/Rosetta_rna_file/VC_riboswitch/create_symmetric_helix_test/NO_mutation_Idealize_BP_alignment_test for proof.

#create_idealize_helix.py -fasta fasta -cutpoint_open 6 -prefect_BP_symmetry_mode True
#create_idealize_helix.py -fasta fasta -prefect_BP_symmetry_mode True


START_argv=copy.deepcopy(argv)

verbose= parse_options( argv, "verbose", "True" )


if(verbose):
	print "----------------------------------------------------------------------------------------------------------------------------------"
	print "----------------------------------------------------------------------------------------------------------------------------------"
	print "Enter: %s " %(list_to_string(START_argv))
	print "----------------------------------------------------------------------------------------------------------------------------------"
	print "----------------------------------------------------------------------------------------------------------------------------------"


HOMEDIR = expanduser('~') 

##########################################################################################
fasta_file = parse_options( argv, "fasta", "" )
input_sequence=parse_options( argv, "sequence", "" ) 

num_input_sequence=0

if(fasta_file!=""): num_input_sequence+=1
if(input_sequence!=""): num_input_sequence+=1

if(num_input_sequence!=1): error_exit_with_message("num_input_sequence=%d!=1" %(num_input_sequence)) 

if(fasta_file!=""):
	if(exists(fasta_file)==False): error_exit_with_message("fasta_file (%s) doesn't exist!!" %(fasta_file)) 
	sequence = open( fasta_file  ).readlines()[1][:-1]
else:
	sequence	=input_sequence

SEQUENCE=sequence.upper()

if(verbose): print "sequence=%s " %(sequence)
if(verbose): print "SEQUENCE=%s " %(SEQUENCE)

total_res=len(sequence)

if(verbose): print "total_res= %d " %(total_res)


##########################################################################################

strand_1_seq_num= parse_options( argv, "strand_1_seq_num", [0] )
strand_2_seq_num= parse_options( argv, "strand_2_seq_num", [0] )

create_VDW_screener=True
helix_type="BLAH"

if(len(strand_1_seq_num)!=0): #OLD WAY

	helix_length=len(strand_1_seq_num)
	VDW_screener_type= parse_options( argv, "VDW_screener_type", "NONE" )

	if(VDW_screener_type=="NONE"):
		helix_type="LOWER"
		create_VDW_screener=False
	elif(VDW_screener_type=="LOWER" or VDW_screener_type=="UPPER"):
		helix_type=VDW_screener_type
	else:
		error_exit_with_message("Invalid VDW_screener_type (%s) " %(VDW_screener_type) )
	
else: #OLD WAY

	cutpoint_open= parse_options( argv, "cutpoint_open", 0 )
	helix_length= parse_options( argv, "helix_length", 2 ) #choose 2 to be consistent with old code.

	if(cutpoint_open==0): 
		helix_type="LOWER"
		strand_1_seq_num=range(1,helix_length+1)
		strand_2_seq_num=range(total_res-helix_length+1,total_res+1)
	else:
		helix_type="UPPER"
		strand_1_seq_num=range(cutpoint_open-helix_length+1,cutpoint_open+1)
		strand_2_seq_num=range(cutpoint_open+1,cutpoint_open+helix_length+1)
		
if(len(strand_1_seq_num)!=len(strand_2_seq_num)): error_exit_with_message("len(strand_1_seq_num)=%(d)!=(%d)=len(strand_2_seq_num)" %(strand_1_seq_num, strand_2_seq_num) )

if(helix_length>8): error_exit_with_message("helix_length>8!")

##########################################################################################

verbose= parse_options( argv, "verbose", "True" )

helix_filename= parse_options( argv, "helix_filename", "" )
invert_helix= parse_options( argv, "invert_helix", "False" )

if(helix_filename==""): error_exit_with_message('error, user need to specify helix_filename=""!') 

if(helix_filename[-4:]!=".pdb"): helix_filename=helix_filename + ".pdb"

VDW_screener_pdb_name="%s_VDW_rep_screener.pdb" %(helix_filename[0:-4])

if(exists(helix_filename)): error_exit_with_message("helix_filename (%s) already exist! " %(helix_filename) )

if(exists(VDW_screener_pdb_name)):  error_exit_with_message("VDW_screener_pdb_name (%s) already exist! " %(VDW_screener_pdb_name) )


##########################################################################################

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )



#if(prefect_BP_symmetry_mode==True): #Basically want a helix that is symmetric with itself.
#	error_exit_with_message("prefect_BP_symmetry_mode==True, this mode is not CORRECTLY implemented yet. Implementation below is wrong") 
	
ten_AU_BP_pdb = get_PYPATH("database/IDEALIZE_HELIX/10_AU_BP.pdb")

four_BP_5_UG_3_wobble_pdb = get_PYPATH("database/IDEALIZE_HELIX/GENERIC_5_U_cutpoint_G_3_WOBBLE/4BP_5_UG_3_wobble.pdb")

invert_four_BP_5_UG_3_wobble_pdb = get_PYPATH("database/IDEALIZE_HELIX/GENERIC_5_U_cutpoint_G_3_WOBBLE/invert_chain_seq_cutpoint_4_4BP_5_UG_3_wobble.pdb")

###"WARNING: Temporarily using 5'-GU-3' wobble helical stub extract from 2R8S_J5_5A..until have time to generate generic version" 
four_BP_5_GU_3_wobble_pdb = get_PYPATH("database/IDEALIZE_HELIX/2R8S_J5_J5A_3_U_cutpoint_G_5_WOBBLE/2R8S_J5_JA_4BP_5_GU_3_wobble.pdb")

invert_four_BP_5_GU_3_wobble_pdb = get_PYPATH("database/IDEALIZE_HELIX/2R8S_J5_J5A_3_U_cutpoint_G_5_WOBBLE/invert_chain_seq_cutpoint_4_2R8S_J5_JA_4BP_5_GU_3_wobble.pdb")

if(invert_helix==True): error_exit_with_message("This option is temporary disallowed")

storage_folder="create_ideal_helix/"

if( exists(storage_folder) ):
	if(verbose): print "Warning storage_folder (%s) already exist ... removing " %(storage_folder) 
	submit_subprocess("rm -r %s " %(storage_folder) )

submit_subprocess("mkdir %s " %(storage_folder) )


main_folder=os.path.abspath(".")


os.chdir( storage_folder )




helix_template_pdb=""

mutate_res_string=""
element_slice_segment=""
VDW_screener_slice_segment=""
skip_mutate_res=[]


if(helix_type=="LOWER"):

	strand_1_sequence=SEQUENCE[strand_1_seq_num[0]-1:strand_1_seq_num[-1]]
	strand_2_sequence=SEQUENCE[strand_2_seq_num[0]-1:strand_2_seq_num[-1]]

	if(verbose):
		print "strand_1_seq_num=%s strand_2_seq_num=%s " %(list_to_string(strand_1_seq_num), list_to_string(strand_2_seq_num)),
		print " strand_1_sequence=5'-%s-3' strand_2_sequence=5'-%s-3' " %(strand_1_sequence, strand_2_sequence)

	#not invert_helix	
	if(Is_perfect_WC_helix(strand_1_sequence,strand_2_sequence)):
		helix_template_pdb=ten_AU_BP_pdb

		helix_1_seq_num=range(10-helix_length+1, 10+1)
		helix_2_seq_num=range(10+1, 10+helix_length+1)

		###Feb 03, 2012: Warning, right now the o2star torsion of the nts in the pdb varies from nts to nts so replacing the 
		###helix [9, 10, 11, 12] with the helix [8, 9, 12, 13] might lead to unexpected randomness. Need to be careful here!

	else: 

		if(helix_length!=2): error_exit_with_message("not perfect_Watson_Crick AND helix_length!=2")
	
		helix_1_seq_num=range(4-helix_length+1, 4+1)
		helix_2_seq_num=range(4+1, 4+helix_length+1)

		if(SEQUENCE[strand_1_seq_num[-1]-1]=="U" and SEQUENCE[strand_2_seq_num[0]-1]=="G"):

			skip_mutate_res=[ helix_1_seq_num[-1], helix_2_seq_num[0] ] 
			helix_template_pdb=four_BP_5_UG_3_wobble_pdb

		elif(SEQUENCE[strand_1_seq_num[-1]-1]=="G" and SEQUENCE[strand_2_seq_num[0]-1]=="U"):

			print "WARNING: Temporarily using 5'-GU-3' wobble helical stub extract from 2R8S_J5_5A..until have time to generate generic version" 
			skip_mutate_res=[ helix_1_seq_num[-1], helix_2_seq_num[0] ]
			helix_template_pdb=four_BP_5_GU_3_wobble_pdb

		else:
			error_exit_with_message("A idealized helix cannot be created for strand_1_sequence=%s, strand_2_sequence=%s" %(strand_1_sequence, strand_2_sequence) )

	element_slice_segment=     "%d %d %d %d" %(helix_1_seq_num[0], helix_1_seq_num[-1], helix_2_seq_num[0], helix_2_seq_num[-1])
	VDW_screener_slice_segment="%d %d %d %d" %(helix_1_seq_num[0]-2, helix_1_seq_num[0], helix_2_seq_num[-1], helix_2_seq_num[-1]+2)  

elif(helix_type=="UPPER"):

	strand_1_sequence=SEQUENCE[strand_1_seq_num[0]-1:strand_1_seq_num[-1]]
	strand_2_sequence=SEQUENCE[strand_2_seq_num[0]-1:strand_2_seq_num[-1]]

	if(verbose):
		print "strand_1_seq_num=%s strand_2_seq_num=%s " %(list_to_string(strand_1_seq_num), list_to_string(strand_2_seq_num)),
		print " strand_1_sequence=5'-%s-3' strand_2_sequence=5'-%s-3' " %(strand_1_sequence, strand_2_sequence)

	#not invert_helix	
	if(Is_perfect_WC_helix(strand_1_sequence,strand_2_sequence)):
		helix_template_pdb=ten_AU_BP_pdb

		helix_1_seq_num=range(1, helix_length+1)
		helix_2_seq_num=range((2*10)+1-helix_length, (2*10)+1)

		###Feb 03, 2012: Warning, right now the o2star torsion of the nts in the pdb varies from nts to nts so replacing the 
		###helix [1, 2, 19, 20] with the helix [2, 3, 18, 19] might lead to unexpected randomness. Need to be careful here!

	else:

		if(helix_length!=2): error_exit_with_message("not perfect_Watson_Crick AND helix_length!=2")
	
		helix_1_seq_num=range(1, helix_length+1)
		helix_2_seq_num=range((2*4)+1-helix_length, (2*4)+1)

		if(SEQUENCE[strand_1_seq_num[0]-1]=="G" and SEQUENCE[strand_2_seq_num[-1]-1]=="U"):

			skip_mutate_res=[ helix_1_seq_num[0], helix_2_seq_num[-1] ] 
			helix_template_pdb=invert_four_BP_5_UG_3_wobble_pdb

		elif(SEQUENCE[strand_1_seq_num[0]-1]=="U" and SEQUENCE[strand_2_seq_num[-1]-1]=="G"):

			print "WARNING: Temporarily using 5'-GU-3' wobble helical stub extract from 2R8S_J5_5A..until have time to generate generic version" 
			skip_mutate_res=[ helix_1_seq_num[0], helix_2_seq_num[-1] ] 
			helix_template_pdb=invert_four_BP_5_GU_3_wobble_pdb

		else:
			error_exit_with_message("A idealized helix cannot be created for strand_1_sequence=%s, strand_2_sequence=%s" %(strand_1_sequence, strand_2_sequence) )

	element_slice_segment=     "%d %d %d %d" %(helix_1_seq_num[0], helix_1_seq_num[-1], helix_2_seq_num[0], helix_2_seq_num[-1])
	VDW_screener_slice_segment="%d %d %d %d" %(helix_1_seq_num[-1], helix_1_seq_num[-1]+2, helix_2_seq_num[0]-2, helix_2_seq_num[0])  

else:
	error_exit_with_message("Invalid helix_type (%s) " %(helix_type) )

##############################################################################################################################

mutate_res_string=""

for n in range(len(strand_1_seq_num)): 
	if(helix_1_seq_num[n] in skip_mutate_res): continue
	mutate_res_string+="%s%d-" %(SEQUENCE[strand_1_seq_num[n]-1], helix_1_seq_num[n])

for n in range(len(strand_2_seq_num)): 
	if(helix_2_seq_num[n] in skip_mutate_res): continue
	mutate_res_string+="%s%d-" %(SEQUENCE[strand_2_seq_num[n]-1], helix_2_seq_num[n])

mutate_res_string=mutate_res_string[:-1]

#############################################################################################################################


if(helix_template_pdb==""): error_exit_with_message("helix_template_pdb==\"\"")

if(PATH_exists(helix_template_pdb)==False): error_exit_with_message("helix_template_pdb (%s) doesn't exist!!" %(helix_template_pdb))

if(helix_template_pdb==basename(helix_template_pdb)):
	error_exit_with_message("helix_template_pdb==basename(helix_template_pdb), helix_template_pdb=%s" %(helix_template_pdb))

submit_subprocess("cp %s %s " %(helix_template_pdb, basename(helix_template_pdb)) )

helix_template_pdb=basename(helix_template_pdb)

submit_subprocess("replace_chain_inplace.py  %s A " %(helix_template_pdb) )
#############################################################################################################################

SWA_mutate_command_line="SWA_mutate_residues.py " + helix_template_pdb + " " + mutate_res_string

if(verbose): print "SWA_mutate_command_line %s " %(SWA_mutate_command_line)

submit_subprocess(SWA_mutate_command_line)

mutate_helix_pdb="mutate_" + mutate_res_string + "_" + helix_template_pdb

if(exists(mutate_helix_pdb)==False):  error_exit_with_message("mutate_helix_pdb (%s) doesn't exist!!" %(mutate_helix_pdb))


#######################################################
'''
if(invert_helix==True): #Nov 23, 2010. 
	
	inverted_mutate_helix_pdb="inverted_" + mutate_helix_pdb
	invert_helix_command="invert_chain_sequence.py -cutpoint_open %s -input_pdb_file %s -output_pdb_file %s " %(4, mutate_helix_pdb, inverted_mutate_helix_pdb) 
	#4 is the cutpoint open of the idealized helix and not the actual motif we are trying model (the one corresponding to the fasta file)

	if(verbose): print "invert_helix_command= %s" %(invert_helix_command)
	submit_subprocess(invert_helix_command)

	mutate_helix_pdb=inverted_mutate_helix_pdb
'''
#######################################################

create_element_command="SWA_pdbslice.py %s -segments %s %s" %(mutate_helix_pdb, element_slice_segment, helix_filename)
if(verbose): print "create_element_command= %s" %(create_element_command) 
submit_subprocess(create_element_command)

if(exists(helix_filename)==False): error_exit_with_message("helix_filename=%s doesn't exist!" %(helix_filename) )

submit_subprocess("renumber_pdb_in_place.py %s" %( helix_filename) ) 

#######################################################

if(create_VDW_screener):
	create_VDW_screener_command="SWA_pdbslice.py %s -segments %s %s" %(mutate_helix_pdb, VDW_screener_slice_segment, VDW_screener_pdb_name)
	if(verbose): print "create_VDW_screener_command= %s" %(create_VDW_screener_command) 
	submit_subprocess(create_VDW_screener_command)

	if(exists(helix_filename)==False): error_exit_with_message("VDW_screener_pdb_name=%s doesn't exist!" %(helix_filename) )

	submit_subprocess("renumber_pdb_in_place.py %s" %( VDW_screener_pdb_name) ) 

	#lower_VDW_screener.pdb  3-4 1-16 upper_VDW_screener.pdb  1-6 8-9 '
	VDW_rep_screen_info_string= "-VDW_rep_screen_info %s" %(VDW_screener_pdb_name)
	if(helix_type=="LOWER"):
		VDW_rep_screen_info_string+=" 3-4 %d-%d " %(strand_1_seq_num[0], strand_2_seq_num[-1])
	else:
		VDW_rep_screen_info_string+=" 1-6 %d-%d " %(strand_1_seq_num[-1], strand_2_seq_num[0])


#######################################################

os.chdir( main_folder )

submit_subprocess("cp %s/%s %s " %(storage_folder, helix_filename, helix_filename))

if(create_VDW_screener):
	submit_subprocess(' echo "%s" >> VDW_rep_screen_info.txt' %(VDW_rep_screen_info_string) )
	submit_subprocess("cp %s/%s %s " %(storage_folder, VDW_screener_pdb_name, VDW_screener_pdb_name))





