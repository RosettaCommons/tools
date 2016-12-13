#!/usr/bin/env python

######################################################################
from SWA_DAG_util import *
from SWA_sampling_post_process_util import *
######################################################################

from SWA_dagman_python.misc.SWA_cat_outfiles import *
from SWA_dagman_python.misc.SWA_filter_silent_file import *
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
from os.path import basename, dirname, exists, expanduser
from time import sleep


#############################################################################

python_command=get_python_command(argv)
print_title_text("Enter " + python_command)

#Removeed all score filtering on Jan 12, 2012

indir_prefix = parse_options( argv, "indir_prefix", "" )
if(indir_prefix == ""): error_exit_with_message("need to input indir_prefix!")

condor_submit_file=parse_options( argv, "condor_submit_file", "" ) #FEB 08, 2012
if(condor_submit_file == ""): error_exit_with_message("condor_submit_file == \"\" ") 

if(exists(condor_submit_file)==False): error_exit_with_message("condor_submit_file (%s) doesn't exist!" %(condor_submit_file) )

min_filtered_struct=parse_options(argv, "min_filtered_struct", 8000) #previous default value
print "min_filtered_struct =%d " %(min_filtered_struct)

delete_files=True

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )

################################################################################################

###Update this section on Jan 10, 2012
#filter_outfile = cat_outfile.replace('.out','_filtered.out')  ####BEWARE THIS IS SPECIFIC FORM IS HARD CODED in SWA_rna_build_dagman.py
#cat_outfile = dirname(indir_prefix) + '/' + basename(indir_prefix).lower() + '_sample.out'

cat_outfile = get_sampler_post_process_cat_outfile( indir_prefix )
filter_outfile = get_sampler_post_process_filter_outfile( indir_prefix )

deleting_files_signal_file=get_sampler_post_process_generic_file(indir_prefix, prepend_str="deleting_files_signal_", append_str =".txt")

if(exists(deleting_files_signal_file)):	error_exit_with_message("deleting_files_signal_file (%s) exist!"	%(deleting_files_signal_file))

################################################################################################

if(exists(cat_outfile)): 
	print "Warning cat_outfile (%s) already exist....removing!" %(cat_outfile)
	submit_subprocess("rm %s" %(cat_outfile))

no_non_empty_sampler_silent_file=concatenate_sampler_silent_file(condor_submit_file, cat_outfile)

if(no_non_empty_sampler_silent_file):
	print "EARLY EXIT, no_non_empty_sampler_silent_file! | python_command=%s" %(python_command)
	outfile = open( filter_outfile , 'w')
	outfile.write("empty filtered silent_file since no non-empty sampler silent_file. indir_prefix = %s" %(indir_prefix) )
	outfile.close()

	if(delete_files):	create_generic_done_signal_file(deleting_files_signal_file) 	#The mark the point where the scipt is now not rerunnable!

	delete_sampler_outfiles_and_folders(indir_prefix, delete_files) #possible that empty folders exist, delete them.	

else:

	if(exists(cat_outfile)==False): error_exit_with_message("cat_outfile (%s) doesn't exist!" %(cat_outfile) )

	filter_sampler_cat_outfile(filter_outfile, cat_outfile, min_filtered_struct)

	if(exists(filter_outfile)==False): error_exit_with_message("filter_outfile (%s) doesn't exist!" %(filter_outfile) )

	#############DELETE files and folders that are no longer needed to save DISK-SPACE!###########
	if(delete_files):	create_generic_done_signal_file(deleting_files_signal_file)		#The mark the point where the scipt is now not rerunnable!

	delete_sampler_outfiles_and_folders(indir_prefix, delete_files) #THIS USED TO BE CALLED BEFORE the call to filter_silent_file (Moved down here on Jan 12, 2012)

	if(delete_files): ###Delete the cat_outfile

		if(exists(cat_outfile)==False): error_exit_with_message("cat_outfile (%s) doesn't exist!" %(cat_outfile))

		command = 'rm %s '%(cat_outfile)
		print( command )
		sys.stdout.flush()
		sys.stderr.flush()
		submit_subprocess( command )

####################################################################################################

print_title_text("Exit " + python_command)

