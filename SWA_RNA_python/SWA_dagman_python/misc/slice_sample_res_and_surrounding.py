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
from SWA_dagman_python.parser.SWA_parse_options  import parse_options
######################################################################

#assert( len(argv)>2)
#pose_name = argv[1]
#print pose_name

#pose_name='1s72_RNA_A.pdb'

#June + July, 2011
#slice_sample_res_and_surrounding.py -loop_segment_list 117-121  -pose_name 3d2v_RNA_A.pdb -MODE single_node -expand_radius_list 10.0 30.0 > log.txt
#slice_sample_res_and_surrounding.py -loop_segment_list 286-292 -pose_name 3g78_RNA_A.pdb -MODE single_node -expand_radius_list 10.0 30.0 > log.txt

#slice_sample_res_and_surrounding.py -loop_segment_list 21-25 -pose_name 2r8s_RNA.pdb -MODE single_node -expand_radius_list 10.0 30.0 > log.txt

#slice_sample_res_and_surrounding.py -loop_segment_list 1950-1959  -pose_name 1s72_RNA_A.pdb -MODE single_node -expand_radius_list 10.0 30.0 > output.txt

###NOTE ellipsoid_envelope_mode=False for the jobs BELOW!##########

#slice_sample_res_and_surrounding.py -loop_segment_list 35-40 520-525 1923-1932 1950-1959 2376-2382 -pose_name 1s72_RNA_A.pdb -MODE standard_queue  -expand_radius_list 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 15.0 20.0 25.0 30.0 

#March 07, 2011
#slice_sample_res_and_surrounding.py -loop_segment_list 7-11  -pose_name 2pn4_RNA_A.pdb -MODE single_node  -additional_slice_res 5

#March 07, 2011

#Feb 22, 2011
#slice_sample_res_and_surrounding.py -loop_segment_list 54-57 -pose_name 1nuj_RNA_A.pdb -MODE single_node
#slice_sample_res_and_surrounding.py -loop_segment_list 72-76 -pose_name 2qbz_RNA_A.pdb -MODE single_node
#slice_sample_res_and_surrounding.py -loop_segment_list 72-76 -pose_name 2qbz_RNA_A.pdb -MODE single_node

#slice_sample_res_and_surrounding.py -loop_segment_list 173-178 -pose_name 1u6b_RNA_A.pdb -MODE single_node
#slice_sample_res_and_surrounding.py -loop_segment_list 142-145 -pose_name 2qbz_RNA_A.pdb -MODE single_node -additional_slice_res 147
#Feb 22, 2011

#slice_sample_res_and_surrounding.py -loop_segment_list 117-121  -pose_name 3d2v_RNA_A.pdb -MODE single_node 
#slice_sample_res_and_surrounding.py -loop_segment_list 161-167  -pose_name 3owi_RNA_A.pdb -MODE single_node -additional_slice_res 169
#slice_sample_res_and_surrounding.py -loop_segment_list 61-67  -pose_name combine_2qwy.pdb -MODE standard_queue 
#slice_sample_res_and_surrounding.py -loop_segment_list 203-207 -pose_name 1u6b_RNA_A.pdb -MODE single_node -additional_slice_res 209
#slice_sample_res_and_surrounding.py -loop_segment_list 286-292 -pose_name 3g78_RNA_A.pdb -MODE single_node  #standard_queue

#rsync -avz  sripakpa@biox2.stanford.edu:~/minirosetta/SLICE_3G78_LOOP_286_292/. ~/minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/3G78_LOOP_286_292/
#rsync -avz  sripakpa@biox2.stanford.edu:~/minirosetta/SLICE_2QWY_LOOP_60_66_BIOX_COMBINE/. ~/minirosetta/Rosetta_rna_file/PAPER_BENCHMARK_PRISTINE/2QWY_LOOP_60_66/BIOX_COMBINE/


####################################################################
copy_argv=copy.deepcopy(argv)

####################################################################

no_graphic_string= parse_options( argv, "no_graphic", "" )



if(use_new_src_code()): 
	EXE = get_rosetta_EXE_specified_no_graphic_string("swa_rna_util", no_graphic_string) 
else:
	EXE = get_rosetta_EXE_specified_no_graphic_string("parin_test", no_graphic_string) 

database_folder= get_rosetta_database_folder() 

####################################################################

additional_slice_res=parse_options( argv, "additional_slice_res", [-1])

#start_res_list=[]
#end_res_list=[1932]
MODE=parse_options( argv, "MODE", "") 

if(MODE!="manual_sub" and MODE!="standard_queue" and MODE!="single_node"):
	error_exit_with_message("Invalid MODE(%s)!" %(MODE))



loop_segment_list=parse_options( argv, "loop_segment_list", [""] ) 

ellipsoid_envelope_mode=parse_options( argv, "ellipsoid_envelope_mode", "True" ) 

if(loop_segment_list==[""]): error_exit_with_message("User need to pass in loop_segment_list!")

expand_radius_list= parse_options( argv, "expand_radius_list", [5.0, 30.0] ) #defualt used to be 10.0 for starting and 30.0 for VDW rep (for the ribosomal loops)

pose_name= parse_options( argv, "pose_name", "1s72_RNA_A.pdb" )

pose_name=os.path.abspath(pose_name)

if( exists( pose_name ) ==False):  error_exit_with_message("pose_name (%s) doesn't exist!" %(pose_name) )

###########################consistency check#######################################

start_res_list=[]
end_res_list=[]

for loop_segment in loop_segment_list:

	loop_segment_split=loop_segment.split('-')
	if(len(loop_segment_split)!=2): error_exit_with_message("len(loop_segment_split)!=2" )

	start_res=int(loop_segment_split[0])
	end_res=int(loop_segment_split[1])

	if(start_res<0):  error_exit_with_message("start_res is not a positive integer, start_res=%d" %(start_res) )
	if(end_res<0):  error_exit_with_message("end_res is not positive integer, end_res=%d" %(end_res) )

	if(start_res>=end_res): error_exit_with_message("start_res(%d) > end_res(%d)" %(start_res, end_res))

	start_res_list.append(start_res)
	end_res_list.append(end_res)
##############################################################################


#start_res_list= parse_options( argv, "start_res_list", [0] ) 
#end_res_list= parse_options( argv, "end_res_list", [0] ) 

if(len(start_res_list)!=len(end_res_list)):
	print "start_res_list=", start_res_list
	print "end_res_list=", end_res_list
	error_exit_with_message("len(start_res_list)!=len(end_res_list)" )

for n in range(len(start_res_list)):

	start_res=start_res_list[n]
	end_res=end_res_list[n]

	foldername="Get_surr_%4s_%4s" %(str(start_res).zfill(4), str(end_res).zfill(4))

	if( exists(foldername) ):   error_exit_with_message("foldername (%s) already exist!" %(foldername) )

	submit_subprocess( 'mkdir %s' %(foldername)  )

	os.chdir( foldername)

	
	sample_res=''

	for seq_num in range(start_res, end_res+1):
		sample_res+='%d ' %(seq_num)


	for expand_radius in expand_radius_list:

		expand_radius_subfolder="expand_%d_Angstrom" %(int(expand_radius+0.0001))

		if( exists(expand_radius_subfolder) ):   error_exit_with_message("expand_radius_subfolder (%s) already exist!" (expand_radius_subfolder) )

		submit_subprocess( 'mkdir %s' %(expand_radius_subfolder)  )

		os.chdir( expand_radius_subfolder)

		command=EXE
		if(ellipsoid_envelope_mode):
			command += " -algorithm slice_ellipsoid_envelope"
		else:
			command += " -algorithm slice_sample_res_and_surrounding"
		command += " -database %s" %(database_folder)
		command += " -s " + pose_name
		command += " -sample_res " + sample_res
		command += " -surrounding_radius %s " %(expand_radius)

		if(len(additional_slice_res)!=0): command += " -additional_slice_res %s " %(list_to_string(additional_slice_res))

		#command += " > output.txt"

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

		elif(MODE=="single_node"):	

			command += " > log.out 2> log.err"
			print command
			submit_subprocess( command )

		else:
			error_exit_with_message("Invalid mode=%s "%(MODE))

		os.chdir( "../")

	os.chdir( "../")

print 

#bsub -q IA  -o hello.out -e hello.err -Is SWA_hello_world.py

print "----------------------------------------------------------------------------------------------------------------------------"
print "----------------------%s----------------------" %(list_to_string(copy_argv))
print "slice_sample_res_and_surrounding.py sucessfully RAN! "
print "----------------------------------------------------------------------------------------------------------------------------"

