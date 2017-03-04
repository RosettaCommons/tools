#!/usr/bin/env python

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
from os.path import basename, dirname, exists, expanduser
from time import sleep

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options, replace_arg_value
########################################################################

HOMEDIR = expanduser('~') 

local_dag_file=parse_options( argv, "local_dag_file", "" )
VERBOSE=parse_options( argv, "VERBOSE", "True" )
minimizer_rename_tag=parse_options( argv, "minimizer_rename_tag", "False" )
output_pdb=parse_options( argv, "output_pdb", "True" )
clusterer_silent_file_in=parse_options( argv, "clusterer_silent_file_in", "" )

dag_file = parse_options( argv, "dag_file", "" )
if(dag_file==""): error_exit_with_message("User need pass in dag_file (%s) !!" %(dag_file)) 
if(exists(dag_file)==False): error_exit_with_message("dag_file (%s) doesn't exist!!" %(dag_file)) 

queue_ID = parse_options( argv, "queue_ID", 0 )

folder_layer  = parse_options( argv, "folder_layer", 1 )

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )

JOB_TYPE=""

if(len(dag_file.split("/"))>=3):
	if(dag_file.split("/")[-3]=="CONDOR"):
		if(dag_file.split("/")[-2] in ["CLUSTERER", "DS_COMBINE_CLUSTERER" ,"FILTERER", "REBUILD_BULGE", "SAMPLER" , "VIRT_RIBOSE_SAMPLER"]):
			JOB_TYPE=dag_file.split("/")[-2] + "_"

if(local_dag_file==""):
	local_dag_file="LOCAL_%s%s" %(JOB_TYPE, basename(dag_file))

double_dot_prefix=""

if(folder_layer==0):
	double_dot_prefix=""
elif(folder_layer==1):
	double_dot_prefix="../"
elif(folder_layer==2):
	double_dot_prefix="../../"
else:
	error_exit_with_message("folder_layer (%d) need to equal 1 or 2!" %(folder_layer) ) 


lines = safe_open( dag_file, mode='r' ,Is_master=False ).readlines()

exe= ""
args= ""

for line in lines:
	if len( line ) > 2: 
		cols = string.split( line )
		if cols[0] ==  "executable":
			assert( cols[1] == "=" )
			exe = cols[2]
		elif cols[0] == "arguments":
			assert( cols[1] == "=" )
			args = string.join(cols[2:])



'''
exe_split=exe.split('/')

if(exe_split[1] not in ["home", "Users"]): error_exit_with_message("exe_split[1] not in [\"home\", \"Users\"]")

exe_split[2] is the user name, e.g sripakpa, rhiju

actual_exe=HOMEDIR + "/" + list_to_string(exe_split[3:],'/')[1:]

print "exe        = %s" %(exe)
print "actual_exe = %s" %(actual_exe)
'''

actual_exe=exe

if(is_release_mode()==False): #For local_run on Parin Sripakdeevong Mac laptop!

	######OK this could be dangerous, since the CONDOR file might be out of date + depend on either RELEASE_MODE is True or False###
	#if(actual_exe.split("/")[1]=="home"):
	#if( (actual_exe.split("/")[2]=="sripakpa") or (actual_exe.split("/")[2]=="kwipapat") or (actual_exe.split("/")[2]=="vanlang")):

	if( (actual_exe.split("/")[-3]=="rosetta_source") and (actual_exe.split("/")[-2]=="bin") and (actual_exe.split("/")[-1][-15:]=="linuxgccrelease") ):
		actual_exe=get_rosetta_EXE(actual_exe.split("/")[-1].replace(".linuxgccrelease", "") )

print "actual_exe = %s" %(actual_exe)

if(PATH_exists(actual_exe)==False): error_exit_with_message("actual_exe (%s) doesn't exist!" %(actual_exe))

args_list=args.split()

for n in range(len(args_list)-1, -1, -1):
	curr_arg=args_list[n]
	#if(curr_arg.count("minirosetta_database")==1): args_list[n]="/Users/sripakpa/minirosetta_database" #Comment out on Jan 01, 2012!
	if(curr_arg.count(".pdb")==1								): args_list[n]= double_dot_prefix + curr_arg
	if(curr_arg.count(".out")==1								): args_list[n]= double_dot_prefix + curr_arg
	if(curr_arg.count(".bmrb")==1								): args_list[n]= double_dot_prefix + curr_arg

	if(curr_arg=="-fasta"												): args_list[n+1]= double_dot_prefix + args_list[n+1]
	if(curr_arg=="-chemical_shift_args_file"			): args_list[n+1]= double_dot_prefix + args_list[n+1] #Nov 20, 2011
	if(curr_arg=="-out:file:silent"							): args_list[n+1]= basename(args_list[n+1])
	if(curr_arg=="-output_filename"							): args_list[n+1]= basename(args_list[n+1])
	if(curr_arg=="-filter_output_filename"			): args_list[n+1]= basename(args_list[n+1])
	if(curr_arg=="-pre_process_output_filename"	): args_list[n+1]= double_dot_prefix + args_list[n+1]

	#if(basename(exe)=="SWA_cluster.py"):
	#	if(curr_arg.count("-VERBOSE")==1): args_list[n+1]= "false"


actual_args=""
for arg in args_list:
	actual_args+=" " + arg

if( exe.count(".linuxgccrelease")>0 or exe.count(".macosgccrelease") ):

	if(output_pdb):						actual_args+=" -output_pdb true "

	if(actual_args.count("rna_resample_test")==1 or actual_args.count("rna_sample")==1): 

		if(VERBOSE): 												actual_args=replace_arg_value(actual_args, "VERBOSE", "true", allow_missing=True)
		if(minimizer_rename_tag==False):			actual_args=replace_arg_value(actual_args, "minimizer_rename_tag", "false", allow_missing=True)
		actual_args+= " > rna_sample_output.txt"
	elif(actual_args.count("rna_sample_virtual_ribose")==1):
		actual_args+= " > sample_virtual_ribose_output.txt"
	elif(actual_args.count("filter_combine_long_loop")==1):
		actual_args+= " > filterer_output.txt"
	elif(actual_args.count("clusterer")>0):
		if(clusterer_silent_file_in!=""): 		actual_args=replace_arg_value(actual_args, "in:file:silent", clusterer_silent_file_in, allow_missing=False)
		actual_args+= " > rna_cluster_output.txt"
	elif(actual_args.count("rna_minimize")>0):
		if(VERBOSE): 												actual_args=replace_arg_value(actual_args, "VERBOSE", "true", allow_missing=True)
		if(minimizer_rename_tag==False):			actual_args=replace_arg_value(actual_args, "minimizer_rename_tag", "false", allow_missing=True)
		actual_args+= " > rna_minimize_output.txt"
	else:
		error_exit_with_message("Invalid Rosetta exe (%s)" %(exe) )

elif( exe.count('/SWA_dagman_python/')>0 ):

	if(exe.count("DAG_rebuild_bulge.py")>0): actual_args+= " -location LOCAL "

	actual_args+= " > python_DAG_output.txt"

else:
	error_message("Invalid exe (%s) " %(exe) )


actual_args=actual_args.replace("$(Process)", "%d" %(queue_ID))

#print "actual_exe= ", actual_exe
#print "actual_args= ", actual_args

command= actual_exe + " " + actual_args

print "command= %s" %(command)



submit_subprocess(' echo "\n\n%s\n" > %s' %(command, local_dag_file) )



