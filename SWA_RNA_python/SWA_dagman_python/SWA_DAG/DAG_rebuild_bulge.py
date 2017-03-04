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
from DAG_rebuild_bulge_util import *

######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.utility.DAGMAN_util import parse_chemical_shift_args
from SWA_dagman_python.parser.SWA_parse_options import parse_options, replace_arg_value
################################################################################################################################################################
#DAG_rebuild_bulge.py -silent_file region_FINAL.out -queue_ID 0 -AUTO_REBUILD_RES_AND_CUTPOINT True -main_algorithm Hermann_Phase_Two > log_REBUILD_BULGE.txt

#SWA_rebuild_bulge.py -silent_file cluster_200_score_max_n_struct_500_cat_file.out -start_tag S_0 -AUTO_REBUILD_RES_AND_CUTPOINT True -main_algorithm Hermann_Phase_Two > log_REBUILD_BULGE.txt

START_argv=copy.deepcopy(argv)

print "---------------------------------------------------------------------------------------------------------------------------------" 
print "---------------------------------------------------------------------------------------------------------------------------------" 

start_folder=os.path.abspath(".")

start_time=time.time()

queue_ID= parse_options( argv, "queue_ID", "" ) 
if(queue_ID!=""): 

	try:
		int_queue_ID=int(queue_ID)
	except:
		error_exit_with_message("cannot convert queue_ID (%s) int a integer!" %(queue_ID) )

	start_tag = "S_%s" %(int(queue_ID))

else:

	start_tag= parse_options( argv, "start_tag", "" ) 

if(start_tag==""): error_exit_with_message("start_tag==\"\"")
print "start_tag= %s " %(start_tag) 

start_silent_file = parse_options( argv, "start_silent_file", "" )
if(start_silent_file==""): error_exit_with_message("User need to pass in start_silent_file!")
if(exists( start_silent_file )==False): error_exit_with_message("start_silent_file (%s) doesn't exist!" %(abspath(start_silent_file) ))


main_folder = parse_options( argv, "outfolder", "" ) 
if(main_folder==""): error_exit_with_message("User need to pass in outfolder!")

final_outfile = parse_options( argv, "out:file:silent", "" ) 

common_args_file=parse_options( argv, "common_args_file", "")
if(common_args_file==""): error_exit_with_message("User need to pass in common_args_file!") #Important to check this before parse_sampling_options() line!
if(exists( common_args_file )==False): error_exit_with_message("common_args_file (%s) doesn't exist!" %(abspath(common_args_file) ))
common_args_file=os.path.abspath(common_args_file)


main_algorithm = parse_options( argv, "main_algorithm", "" )
if(main_algorithm!="SWA" and main_algorithm!="FARFAR"): error_exit_with_message("main_algorithm (%s) need to be either \"SWA\" or \"FARFAR\"" %(main_algorithm))

native_pdb =  parse_options( argv, "native_pdb", "" )
if(native_pdb!=""): native_pdb=os.path.abspath(native_pdb) #Important to check this before parse_sampling_options() line!

location= parse_options( argv, "location", "BIOX" )
if(location!="BIOX" and location!="LOCAL"): error_exit_with_message("Invalid location (%s) " %(location))

user_fixed_res_list = parse_options( argv, "fixed_res", [0] )

copy_start_score_line=parse_options( argv, "copy_start_score_line", "True")

BMRB_chemical_shift_file=parse_options( argv,  "BMRB_chemical_shift_file", "")

chemical_shift_args_file = parse_options( argv, "chemical_shift_args_file", "" )

sampling_options=parse_sampling_options(argv, native_pdb)

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )

###########################################################################################################

start_silent_file=os.path.abspath(start_silent_file)

if(final_outfile==""): error_exit_with_message("User need to pass in final_outfile!")
if(exists(final_outfile)): error_exit_with_message("final_outfile (%s) already exist!" %(final_outfile))
final_outfile=os.path.abspath(final_outfile)	


if(BMRB_chemical_shift_file!=""):

	if(exists(BMRB_chemical_shift_file)==False): error_exit_with_message("BMRB_chemical_shift_file (%s) doesn't exist!" %(BMRB_chemical_shift_file) )
	BMRB_chemical_shift_file=os.path.abspath(BMRB_chemical_shift_file)

	if(exists(chemical_shift_args_file)==False): error_exit_with_message("chemical_shift_args_file (%s) doesn't exist!" %(chemical_shift_args_file) )
	chemical_shift_args_file=os.path.abspath(chemical_shift_args_file)

#######################################################

main_folder=os.path.abspath(main_folder)

if(exists(main_folder)):
	print "warning...main_folder:%s already exist...removing it...! " %(main_folder) 
	submit_subprocess("rm -r %s " %(main_folder) )

submit_subprocess("mkdir -p %s" %(main_folder))

os.chdir( main_folder )
#######################################################
#Note that final_outfile should not be in the working_folder!

working_folder="WORKING_FOLDER/"

working_folder=os.path.abspath(working_folder)

if(exists(working_folder)): error_exit_with_message("working_folder (%s) already exist!" %(working_folder) )

submit_subprocess("mkdir -p %s" %(working_folder))

os.chdir( working_folder )
#######################################################

(ANNOTATED_SEQUENCE, start_score_line, start_column_line)= extract_start_tag_info(start_silent_file, start_tag)

(seq_list, virtual_res_list, cutpoint_lower_list)=parse_annotated_sequence(ANNOTATED_SEQUENCE)

total_res=len(seq_list)

#########COMMON_ARGS OPTIONS (More in run_rebuild_bulge_sampling function)##################################################
common_args=safe_readlines( common_args_file )[0].split()

cutpoint_open_list= parse_options( common_args, "cutpoint_open", [0]) #if cutpoint open!=0, then double strand
if(len(cutpoint_open_list)>1): error_exit_with_message("len(cutpoint_open_list)>1, this option have not been tested yet.")

global_sample_res_list=parse_options( common_args, "global_sample_res_list", [0])
if(len(global_sample_res_list)==0): error_exit_with_message("global_sample_res_list option doesn't exist in common_args: %s!" %(list_to_string(common_args)) )

calc_chemical_shift_res_list=add_boundary_seq_num_to_list(global_sample_res_list, total_res, boundary_size=1) 

common_args=""

sampling_options+=extract_sampling_options_from_common_args(common_args_file, virtual_res_list, location, start_folder)
################################################################################################################################################################

print "total_res= %s "  %(total_res)
print "seq_list= " , seq_list
print "virtual_res_list= " , virtual_res_list
print "cutpoint_open_list=  " , cutpoint_open_list
print "cutpoint_lower_list= " , cutpoint_lower_list

rebuild_seq_num_list=virtual_res_list

#cutpoint_closed_list=figure_out_cutpoint_closed_list_old(main_algorithm, rebuild_seq_num_list, cutpoint_open_list, cutpoint_lower_list, total_res)

cutpoint_closed_list=figure_out_cutpoint_closed_list(main_algorithm, rebuild_seq_num_list, start_silent_file, start_tag, total_res)

#empty these list since they are no longer used after this point!
virtual_res_list=[] 
cutpoint_open_list=[]
cutpoint_lower_list=[]

###############Create fasta file############

fasta_file='auto_fasta'

seq_string=""
for char in seq_list: seq_string += char

submit_subprocess('echo ">full_pose" > %s' %(fasta_file) )
submit_subprocess("echo %s >> %s " %(seq_string, fasta_file) )

#########################Submit sampling job##################################

rebuild_struct_exist=True

input_silent_file=start_silent_file
input_tag=start_tag

for bulge_ID in range(len(rebuild_seq_num_list)):

	cutpoint_closed=cutpoint_closed_list[bulge_ID]
	rebuild_seq_num=rebuild_seq_num_list[bulge_ID]


	rebuild_silent_file=run_rebuild_bulge_sampling(input_silent_file, input_tag, rebuild_seq_num, cutpoint_closed, user_fixed_res_list, \
														  									 total_res, fasta_file, location, sampling_options)


	if(exists(rebuild_silent_file)==False): error_exit_with_message("rebuild_silent_file (%s) doesn't exist!" %(rebuild_silent_file))

	silent_data=safe_open(rebuild_silent_file, mode='r', Is_master=False)
	
	first_silent_line=silent_data.readline();

	silent_data.close()

	if( first_silent_line=="StepWiseRNA_Minimizer:: num_pose_outputted==0, empty silent_file!\n" or
	    first_silent_line=="StepWiseMinimizer:: num_pose_outputted==0, empty silent_file!\n" ):

		print "WARNING rebuild_silent_file (%s) contains NO structure ... setting rebuild_struct_exist to False" %(rebuild_silent_file)
		rebuild_struct_exist=False
		break

	else:
		assert_is_valid_non_empty_silent_file(rebuild_silent_file)		

	#if(exists(rebuild_silent_file)==False):
	#	print "WARNING rebuild_silent_file (%s) doesn't exist! ... setting output_struct_exist to False" %(rebuild_silent_file)
	#	rebuild_struct_exist=False
	#	break

	sorted_silent_file=sort_silent_file_by_score(rebuild_silent_file)

	##############################################################################################################################

	if(BMRB_chemical_shift_file!=""):

		BMRB_silent_file="WITH_CHEM_SHIFT_"+ basename(sorted_silent_file)

		if(exists(BMRB_silent_file)): 
			error_exit_with_message("BMRB_silent_file (%s) already exist!" %(BMRB_silent_file) )

		add_chemical_shift_data_command="add_chemical_shift_score_to_silent_file.py -input_BMRB %s " %(BMRB_chemical_shift_file)
		add_chemical_shift_data_command+=" -residue_list %s " %(list_to_string(calc_chemical_shift_res_list))
		add_chemical_shift_data_command+=" -input_silent_file %s -output_silent_file %s " %(sorted_silent_file, BMRB_silent_file)

		add_chemical_shift_data_command+= parse_chemical_shift_args(chemical_shift_args_file)


		print "add_chemical_shift_data_command=%s " %(add_chemical_shift_data_command)

		submit_subprocess(add_chemical_shift_data_command)

		if(exists(BMRB_silent_file)==False): error_exit_with_message("BMRB_silent_file (%s) doesn't exist!" %(BMRB_silent_file) )

		sorted_silent_file=sort_silent_file_by_score(BMRB_silent_file)

	################################################################################################################################

	#update input silent_file and input tag:
	input_tag="S_000000"
	input_silent_file=sorted_silent_file

############################################################################################################################################################

if(exists(final_outfile)): error_exit_with_message("final_outfile (%s) already exist!" %(final_outfile))

if(rebuild_struct_exist==False):	

	print "DAG_builge_bulge.py: Unsucessfully rebuild bulge(s) for %s!!" %(input_tag)

	outfile = open( final_outfile , 'w')
	outfile.write("DAG_builge_bulge.py: Unsucessfully rebuild bulge(s)!!\n" )
	outfile.close()


else: #OK THIS ACCOUNTS FOR THE CASE WHERE THERE IS NO BULGE AS WELL!

	#####Could seperate the case where len(rebuild_seq_num_list)==0#######
	assemble_command = RNA_SWA_MAIN_EXE + ' -algorithm post_rebuild_bulge_assembly '  
	assemble_command+= " -database %s " %(DATABASE_FOLDER)

	if(native_pdb!=""): assemble_command+= " -native " + native_pdb

	assemble_command+= " -start_tag %s " %(start_tag) 

	assemble_command+= " -start_silent %s " %(start_silent_file)

	#################################################################

	assemble_command+= " -out:file:silent %s " %(final_outfile) 

	####################Only 1 struct for now!!#####################

	assemble_command+= " -in:file:tags %s " %(input_tag) 

	assemble_command+= " -in:file:silent %s " %(input_silent_file)

 #################################################################
	common_args_string=safe_readlines( common_args_file )[0]

	common_args_string=replace_arg_value(	common_args_string, "fasta", fasta_file,  allow_missing=False)

	assemble_command+= " %s " %(common_args_string)

	if(location=="LOCAL"):
		assemble_command+= " -VERBOSE true"
		assemble_command+= " -output_pdb true "
	else:
		assemble_command+= " -VERBOSE false"
		assemble_command+= " -output_pdb false "

	assemble_command+= ' > rosetta_assemble_outfile.txt 2> rosetta_assemble_errfile.txt' 

	submit_subprocess(assemble_command)

	if(exists(final_outfile)==False): error_exit_with_message("final_outfile (%s) doesn't exist!" %(final_outfile) )


	##############################################################################################

	##############################################################################################

	if(copy_start_score_line):

		num_final_score_line=0
		num_final_column_line=0
		final_score_line=""

		FINAL_LINES=safe_readlines( final_outfile )

		submit_subprocess("rm %s " %(final_outfile) )

		OUTFILE=open(final_outfile, 'w') #Replace the existing find!

		for line in FINAL_LINES:

			if( line.count('SCORE:') != 0): 
						
				if( line.count('description') != 0 ): #Description Column

					final_column_line=line
					num_final_column_line+=1

					print "-----------------------------------------Replacing final_column_line:-------------------------------------------------"
					print final_column_line[:-1]
					print "-----------------------------------------WITH start_column_line:------------------------------------------------------"
					print start_column_line[:-1]
					print "----------------------------------------------------------------------------------------------------------------------"

					OUTFILE.write(start_column_line)

				else:

					if( (line.split()[-1])!=start_tag ): error_exit_with_message("(line.split()[start_desc_column])!=start_tag")

					final_score_line=line
					num_final_score_line+=1

					print "-----------------------------------------Replacing final_score_line:-------------------------------------------------"
					print final_score_line[:-1]
					print "-----------------------------------------WITH start_score_line:------------------------------------------------------"
					print start_score_line[:-1]
					print "---------------------------------------------------------------------------------------------------------------------"

					OUTFILE.write(start_score_line)

			else:
				OUTFILE.write(line)


	if(num_final_column_line!=1): error_exit_with_message("num_final_column_line!=1 for final_outfile (%s) !" %(final_outfile))

	if(num_final_score_line!=1): error_exit_with_message("num_final_score_line!=1 for start_tag (%s) !" %(start_tag))

	


os.chdir( main_folder )
if(exists(working_folder)==False): error_exit_with_message("working_folder (%s) doesn't exist!" %(working_folder) )
if(location!="LOCAL"): submit_subprocess("rm -r %s " %(working_folder))


print "---------------------------------------------------------------------------------------------------------------------------------" 
print "---------------------------------------------------------------------------------------------------------------------------------" 

print "-------------------------------------------------------------------------------------------"
print "Successfully RAN: %s" %(list_to_string(START_argv))
print "-------------------------------------------------------------------------------------------"

