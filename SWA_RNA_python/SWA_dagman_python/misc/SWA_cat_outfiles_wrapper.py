#!/usr/bin/env python


from sys import argv,exit,stdout,stderr
import sys
from os.path import exists
import string
import time
from glob import glob
######################################################################

from  SWA_dagman_python.misc.SWA_cat_outfiles import *
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################

#SWA_cat_outfiles_wrapper.py -add_file_num_to_tag False -outfile WITH_SHIFT_STATS_UNIQUE_H5_Feb_10_SWA_Dec_01_1JJ2.out -infile_list  replace_scoreline_WITH_SHIFT_STATS_UNIQUE_H5_Feb_10_SWA.out  replace_scoreline_WITH_SHIFT_STATS_UNIQUE_H5_Dec_01_1JJ2.out   


#SWA_cat_outfiles_wrapper.py -add_file_num_to_tag False -outfile WITH_SHIFT_STATS_combine.out -infile_list replace_scoreline_WITH_SHIFT_STATS_Jan_13_1JJ2_rebuild_bulge_clustered_FINAL_filtered_energy.out replace_scoreline_WITH_SHIFT_STATS_Jan_12_RNA_09_rebuild_bulge_clustered_FINAL_filtered_energy.out  replace_scoreline_WITH_SHIFT_STATS_Oct_22_SWA_region_FINAL.out


#SWA_cat_outfiles_wrapper.py -add_file_num_to_tag False -outfile WITH_SHIFT_STATS_combine.out -infile_list replace_scoreline_WITH_SHIFT_STATS_Dec_20_RNA09_rebuild_bulge_clustered_FINAL_filtered_energy.out replace_scoreline_WITH_SHIFT_STATS_Dec_20_SWA_rebuild_bulge_region_FINAL.out

#SWA_cat_outfiles_wrapper.py -add_file_num_to_tag False -outfile combine.out -infile_list replace_scoreline_rebuild_bulge_region_FINAL.out replace_scoreline_rebuild_bulge_clustered_FINAL_filtered_energy.out replace_scoreline_Dec_02_RNA_09_FARFAR_rebuild_bulge_clustered_FINAL_filtered_energy.out

#SWA_cat_outfiles_wrapper.py -outfile RMSD_test.out -infile_pattern DAG_ID_^_filtered_RMSD.out

#SWA_cat_outfiles_wrapper.py -outfile ENERGY_test.out -infile_pattern DAG_ID_^_filtered_energy.out


outfile= parse_options( argv, "outfile", "" )
infile_list=parse_options(argv, "infile_list", [""])
infile_pattern=parse_options(argv, "infile_pattern", "" )

add_file_num_to_tag=parse_options(argv, "add_file_num_to_tag", "True")

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )

if(outfile==""): error_exit_with_message('user need to specify outfile') 

if(infile_list==[""] and infile_pattern==""): error_exit_with_message('user need to specify either infile_list or infile_pattern !') 

if(infile_list!=[""] and infile_pattern!=""): error_exit_with_message('both infile_list and infile_pattern specified!') 


infile_pattern=infile_pattern.replace("^", "*")

if(infile_pattern!=""):
	print "infile_pattern= " , infile_pattern
	infile_list = glob( infile_pattern )	

print "infile_list= ", infile_list

for infile in infile_list:
	if(exists(infile)==False): error_exit_with_message("infile (%s) doesn't exist!" %(infile))


concatenate_outfiles(infile_list, outfile, add_file_num_to_tag)





