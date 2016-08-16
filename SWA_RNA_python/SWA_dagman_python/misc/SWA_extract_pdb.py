#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
import os
from time import sleep
######################################################################

from SWA_dagman_python.parser.SWA_parse_options import parse_options, parse_seq_num_list_option
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.utility.Rosetta_to_standard_PDB_functions import *
from SWA_dagman_python.utility.PDB_operations import renumber_pdb_in_place_func
######################################################################

#sys.path.append("desired directory")

#SWA_extract_pdb.py -silent_file select_models_suite_clustering_cluster_200_WITH_SHIFT_STATS_combine.out -convert_to_STANDARD_PDB True

#SWA_extract_pdb.py -silent_file select_models_suite_clustering_cluster_200_WITH_SHIFT_STATS_combine.out -convert_to_STANDARD_PDB True -slice_res_list 2-9


##############################################################################################################################

no_graphic_string= parse_options( argv, "no_graphic", "" )

HOMEDIR = expanduser('~') 

print "HOMEDIR= %s" %(HOMEDIR)

EXE=get_rosetta_EXE_specified_no_graphic_string("rna_extract", no_graphic_string)

database_folder= get_rosetta_database_folder()

##############################################################################################################################


print
print


silent_file = parse_options( argv, "silent_file", "" )

if(silent_file==""): error_exit_with_message("User need to pass in silent_file!")

if(exists( silent_file )==False): error_exit_with_message("silent_file %s doesn't exist!" %(abspath(silent_file) ))

#remove_variant_types = parse_options( argv, "remove_variant_types", "True" )

remove_variant_cutpoint_atoms=parse_options( argv, "remove_variant_cutpoint_atoms", "True" )

move_silent_file_into_folder = parse_options( argv, "move_silent_file_into_folder", "True" )
output_foldername = parse_options( argv, "output_foldername", "")
convert_to_STANDARD_PDB = parse_options( argv, "convert_to_STANDARD_PDB", "False" )
slice_segment_string_list = get_segment_string_list( parse_seq_num_list_option( argv, "slice_res_list") )

print
print


desc_column=get_description_column(silent_file, verbose=True)

################################################################################################################################################################


tags_list = parse_options( argv, "tag", [""] )

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )

remove_old_folder=False

if(tags_list==[""]): #extract every structure in the silent_file.

	tags_list=[] #Nov 19, 2011 | Remove the empty string!

	remove_old_folder=True

	for line in open(silent_file):

		Is_score_line_method_1=(line.count('SCORE:') == 1)

		Is_score_line_method_2=(line[0:6]=='SCORE:')

		if(Is_score_line_method_1!=Is_score_line_method_2): error_exit_with_message("Is_score_line_method_1=(%s)!=(%s)=Is_score_line_method_2 for line=%s" %(Is_score_line_method_1, Is_score_line_method_2, line))
		
		if(Is_score_line_method_1==False): continue
			
		if(line.find('description') != -1): continue #Line is not a column_name line

		tag=line.split()[-1] #Oct 19, 2010....assume that description is located at last column..change to this because some silent_file have structs with different score column length

		if( desc_column!=(len(line.split())-1) ): error_exit_with_message("WARNING desc_column=(%s)!=(%s)=(len(line.split())-1), for tag (%s) " %(desc_column, len(line.split())-1, tag) )

		tags_list.append(tag)

if(len(tags_list)==0): error_exit_with_message("len(tags_list)==0")
################################################################################################################################################################

#output_foldername="pose_" + basename(silent_file)
base_directory=os.getcwd()
print base_directory


if(move_silent_file_into_folder):
	if(output_foldername==""):
		output_foldername="pose_" + silent_file.replace('/','_')

	if(remove_old_folder):  
		if(exists(	output_foldername)==True ): submit_subprocess( 'rm -r %s' % (output_foldername) )

	if(exists(	output_foldername)==False): submit_subprocess( 'mkdir %s' % (output_foldername) )

	os.chdir( base_directory + '/' + output_foldername )

	if(dirname(silent_file)==""): silent_file='../' + silent_file

#####Aug 11th, 2011. Check that the output pdb_file doesn't exist YET!#################################
for tag in tags_list:
	pdb_filename="%s.pdb" %(tag)
	if(exists(pdb_filename)): error_exit_with_message("pdb_filename (%s) already exist!" %(pdb_filename) )
####################################################################################################
max_structs_per_extraction=10000 #Basically, the code clashes it there are too many structures.

num_extraction=(len(tags_list)/max_structs_per_extraction)

if( (len(tags_list) % max_structs_per_extraction ) !=0 ): num_extraction+=1

print "num_extraction=%s | len(tags_list)=%s | max_structs_per_extraction=%s " %(num_extraction, max_structs_per_extraction, max_structs_per_extraction)

for extraction_ID in range(num_extraction): #Implement multiple extractions on Feb 06, 2012

	min_tag_num=0   +(extraction_ID*max_structs_per_extraction)
	max_tag_num=9999+(extraction_ID*max_structs_per_extraction)

	print "extraction_ID=%s | min_tag_num=%s | max_tag_num=%s" %(extraction_ID, min_tag_num, max_tag_num)

	command=EXE
	command += " -database %s " %(database_folder)

	command += " -tags " + list_to_string(tags_list[min_tag_num:max_tag_num+1])
	command += " -in:file:silent " + silent_file 
	command += " -in:file:silent_struct_type  binary_rna"

	###OLD SRC###
	##./option.cc.gen.hh:option.add( core::options::OptionKeys::out::file::output_virtual, "Output virtual atoms in output of PDB" ).def(false);
	###NEW SRC###
	##./option.cc.gen.hh:option.add( basic::options::OptionKeys::out::file::output_virtual, "Output virtual atoms in output of PDB" ).def(false);

	###########################################
	command += " -output_virtual true "
	if(remove_variant_cutpoint_atoms==True):
		command += " -remove_variant_cutpoint_atoms true " 	 
	###########################################

	###########################################
	'''
	if(remove_variant_types==True):
		###Basically, the goal to remove all variant_types that virtualize atoms.
		###This accomplishes three things:
		###1. The "extra" virtual atoms (OVU1, OVU2, OVL1) are removed.
		###2. The "real" virtual atoms (atoms in virtual_res or virtual_ribose) are unvirtualized and outputted.
		###3. No virtual atoms left in the pose.

		command += " -remove_variant_types true " 
	else:
		####In this condition, all atoms are outputted as is.

		command += " -remove_variant_types false " 
		command += " -output_virtual true " 
	'''
	###########################################

	command += ' > ' + 'extraction_output_%d.txt' %(extraction_ID)
	print command 
	submit_subprocess( command )

###################################################################

for tag in tags_list:
	pdb_filename="%s.pdb" %(tag)
	if(exists(pdb_filename)==False): error_exit_with_message("pdb_filename (%s) doesn't exist!" %(pdb_filename) )

################################################################################################################################################

if(len(slice_segment_string_list)!=0):

	BEFORE_SLICE_PDB_folder="BEFORE_SLICE_PDB/"
	submit_subprocess("mkdir -p %s" %(BEFORE_SLICE_PDB_folder))

	for tag in tags_list:

		pdb_filename="%s.pdb" %(tag)
		if(exists(pdb_filename)==False): error_exit_with_message("pdb_filename (%s) doesn't already exist!" %(pdb_filename) )

		submit_subprocess("mv %s %s/%s" %(pdb_filename, BEFORE_SLICE_PDB_folder, pdb_filename))

		before_slice_pdb_filename="%s/%s" %(BEFORE_SLICE_PDB_folder, pdb_filename)
		if(exists(before_slice_pdb_filename)==False): error_exit_with_message("before_slice_pdb_filename (%s) doesn't already exist!" %(pdb_filename) )

		after_slice_pdb_filename="AFTER_SLICE_%s.pdb" %(tag)

		if(exists(after_slice_pdb_filename)): error_exit_with_message("after_slice_pdb_filename (%s) already exist!" %(after_slice_pdb_filename) )

		slice_segment_string=""
		for slice_segment in slice_segment_string_list:
			if(len(slice_segment.split("-"))!=2): error_exit_with_message("len(slice_segment.split(\"-\"))!=2")
			slice_segment_string += "%d %d" %(int(slice_segment.split("-")[0]), int(slice_segment.split("-")[1])  ) 

		#print "slice_segment_string=%s" %(slice_segment_string) 

		submit_subprocess("SWA_pdbslice.py -segments %s %s %s >> SLICE_LOG.txt" %(slice_segment_string, before_slice_pdb_filename , after_slice_pdb_filename) )

		if(exists(after_slice_pdb_filename)==False): error_exit_with_message("after_slice_pdb_filename (%s) doesn't exist!" %(after_slice_pdb_filename) )

		if(exists(pdb_filename)): error_exit_with_message("pdb_filename (%s) already exist!" %(pdb_filename) )

		renumber_pdb_in_place_func(after_slice_pdb_filename)	

		submit_subprocess( "mv %s %s " %(after_slice_pdb_filename, pdb_filename) )

####################################################################################################################################################

if(convert_to_STANDARD_PDB):

	ROSETTA_PDB_folder="ROSETTA_PDB/"
	submit_subprocess("mkdir -p %s" %(ROSETTA_PDB_folder))

	for tag in tags_list:

		rosetta_pdb_filename="%s.pdb" %(tag)
		if(exists(rosetta_pdb_filename)==False): error_exit_with_message("rosetta_pdb_filename (%s) doesn't already exist!" %(rosetta_pdb_filename) )

		submit_subprocess("mv %s %s/%s" %(rosetta_pdb_filename, ROSETTA_PDB_folder, rosetta_pdb_filename))
		rosetta_pdb_filename="%s/%s" %(ROSETTA_PDB_folder, rosetta_pdb_filename)
		if(exists(rosetta_pdb_filename)==False): error_exit_with_message("rosetta_pdb_filename (%s) doesn't already exist!" %(rosetta_pdb_filename) )

		standard_pdb_filename="%s.pdb" %(tag)
		if(exists(standard_pdb_filename)): error_exit_with_message("standard_pdb_filename (%s) already exist!" %(standard_pdb_filename) )

		Rosetta_to_standard_PDB_func(input_pdb_file=rosetta_pdb_filename, remove_hydrogen=False, output_pdb_file=standard_pdb_filename, VERBOSE=True)

		if(exists(standard_pdb_filename)==False): error_exit_with_message("standard_pdb_filename (%s) doesn't exist!" %(standard_pdb_filename) )


if(move_silent_file_into_folder): os.chdir( base_directory )
