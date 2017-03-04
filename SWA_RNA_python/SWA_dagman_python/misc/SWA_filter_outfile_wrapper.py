#!/usr/bin/env python

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
from os.path import basename, dirname, exists, expanduser
from time import sleep
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.misc.SWA_cat_outfiles import *
from SWA_dagman_python.misc.SWA_filter_silent_file import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################

#SWA_filter_outfile_wrapper.py -infile cluster_0_FINAL_filtered_energy.out -max_n_struct 20000 -scorecol_name score 

#SWA_filter_outfile_wrapper.py -infile recal_rmsd_WITH_SHIFT_STATS_rebuild_bulge_clustered_FINAL_filtered_energy.out -abs_score_cut 2.00  -scorecol_name NEW_Full_L_rmsd -reverse False

#SWA_filter_outfile_wrapper.py -infile combine.out -abs_score_cut 1.50  -scorecol_name NEW_O_loop_rmsd

#SWA_filter_outfile_wrapper.py -infile WITH_SHIFT_STATS_combine.out -abs_score_cut 0.233  -scorecol_name shift_RMSD

#SWA_filter_outfile_wrapper.py -infile WITH_SHIFT_STATS_combine.out  -max_n_struct 1000  -scorecol_name shift_RMSD


#SWA_filter_outfile_wrapper.py -infile start_from_region_0_4_and_7_0_sample_filtered.out -max_n_struct 25000 -scorecol_name score 

#SWA_filter_outfile_wrapper.py -infile start_from_region_0_5_and_8_0_sample_filtered.out -max_n_struct 25000 -scorecol_name score 

#SWA_filter_outfile_wrapper.py -infile start_from_region_7_3_sample_filtered.out -max_n_struct 25000 -scorecol_name score 

#SWA_filter_outfile_wrapper.py -infile start_from_region_7_4_sample_filtered.out -max_n_struct 25000 -scorecol_name score 

#SWA_filter_outfile_wrapper.py -infile start_from_region_8_5_sample_filtered.out -max_n_struct 25000 -scorecol_name score 

#SWA_filter_outfile_wrapper.py -infile start_from_region_9_5_sample_filtered.out -max_n_struct 25000 -scorecol_name score 

##########################

#SWA_filter_outfile_wrapper.py -infile start_from_region_0_6_and_9_0_sample_filtered.out -max_n_struct 25000 -scorecol_name score 

#SWA_filter_outfile_wrapper.py -infile start_from_region_8_5_sample_filtered.out -max_n_struct 25000 -scorecol_name score 

#SWA_filter_outfile_wrapper.py -infile start_from_region_0_5_and_8_0_sample_filtered.out -max_n_struct 25000 -scorecol_name score 

#SWA_filter_outfile_wrapper.py -infile start_from_region_8_4_sample_filtered.out -max_n_struct 25000 -scorecol_name score 

#SWA_filter_outfile_wrapper.py -infile start_from_region_9_6_sample_filtered.out -max_n_struct 25000 -scorecol_name score 


python_command=get_python_command(argv)
print_title_text("Enter " + python_command)

#######################################################################################
score_diff_cut = parse_options( argv, "score_diff_cut", 0.0 )

max_n_struct= parse_options( argv, "max_n_struct", 0 )

abs_score_cut= parse_options( argv, "abs_score_cut", "FALSE" )

num_filter_specified=0

if(score_diff_cut>0.00001): 
	print "User specified score_diff_cut (%s)" %(score_diff_cut)
	num_filter_specified+=1

if(max_n_struct!=0):
	print "User specified max_n_struct (%s)" %(max_n_struct)  
	num_filter_specified+=1

if(abs_score_cut!="FALSE"):
	print "User specified abs_score_cut (%s)" %(abs_score_cut)  
	num_filter_specified+=1

if(num_filter_specified!=1): error_exit_with_message( "num_filter_specified (%s) != 1" %(num_filter_specified) )

#######################################################################################

filter_outfile= parse_options( argv, "filter_outfile", "" )

remove_SCORE_file = parse_options( argv, "remove_SCORE_file", "False" )

scorecol_name = parse_options( argv, "scorecol_name", "score" )

silent_file = parse_options( argv, "infile", "" )

if(silent_file==""):  error_exit_with_message('silent_file=""!')

if(exists(silent_file)==False): error_exit_with_message("silent_file (%s) doesn't exist!" %(silent_file))

reverse = parse_options( argv, "reverse", "True" )

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )
#######################################################################################

if(score_diff_cut>0.00001):

	print "scorecol_name=%s, score_diff_cut=%s " %(scorecol_name, score_diff_cut)

	if(filter_outfile==""): 
		filter_outfile="%s_score_diff_cut_%s_%s" %(scorecol_name , int(score_diff_cut*100), silent_file) 
		if(reverse==False): filter_outfile="REVERSE_" + filter_outfile

	if(exists(filter_outfile)):
		print "filter_outfile (%s) already exist ... removing" %(filter_outfile)
		submit_subprocess("rm %s " %(filter_outfile) )

	(actual_num_structures, actual_score_gap)=filter_silent_file(silent_file=silent_file, outfile_name=filter_outfile, scorecol_name=scorecol_name, REVERSE=reverse, max_score_gap=score_diff_cut, remove_SCORE_file=remove_SCORE_file)

elif(max_n_struct!=0):

	print "scorecol_name=%s, max_n_struct=%s " %(scorecol_name, max_n_struct)

	if(filter_outfile==""): 
		filter_outfile="%s_max_n_struct_%s_%s" %(scorecol_name , max_n_struct, silent_file) 
		if(reverse==False): filter_outfile="REVERSE_" + filter_outfile

	if(exists(filter_outfile)):
		print "filter_outfile (%s) already exist ... removing" %(filter_outfile)
		submit_subprocess("rm %s " %(filter_outfile) )

	(actual_num_structures, actual_score_gap)=filter_silent_file(silent_file=silent_file, outfile_name=filter_outfile, scorecol_name=scorecol_name, REVERSE=reverse, max_n_struct=max_n_struct, remove_SCORE_file=remove_SCORE_file)

elif(abs_score_cut!="FALSE"):

	print "scorecol_name=%s, abs_score_cut=%s " %(scorecol_name, abs_score_cut)

	if(filter_outfile==""): 
		if(float(abs_score_cut)>1.0):
			filter_outfile="%s_abs_score_cut_%s_%s" %(scorecol_name , int((float(abs_score_cut)+0.0000001)*100), silent_file) 
		else:
			filter_outfile="%s_abs_score_cut_0%s_%s" %(scorecol_name , int((float(abs_score_cut)+0.0000001)*1000), silent_file) 

		if(reverse==False): filter_outfile="REVERSE_" + filter_outfile

	if(exists(filter_outfile)):
		print "filter_outfile (%s) already exist ... removing" %(filter_outfile)
		submit_subprocess("rm %s " %(filter_outfile) )

	(actual_num_structures, actual_score_gap)=filter_silent_file(silent_file=silent_file, outfile_name=filter_outfile, scorecol_name=scorecol_name, REVERSE=reverse, abs_score_cut=abs_score_cut, remove_SCORE_file=remove_SCORE_file)

else:
	error_exit_with_message('SHOULD NOT REACH THIS POINT OF THE SCRIPT!')

print "actual_num_structures=%d, actual_score_gap=%f " %(actual_num_structures, actual_score_gap)

sys.stdout.flush()
sys.stderr.flush()


