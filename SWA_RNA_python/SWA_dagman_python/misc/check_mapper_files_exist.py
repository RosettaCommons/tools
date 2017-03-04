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

from SWA_dagman_python.misc.SWA_cat_outfiles import concatenate_outfiles
######################################################################


from SWA_dagman_python.parser.SWA_parse_options import parse_options
from SWA_dagman_python.utility.SWA_util import *
######################################################################

###############################12:43 AM Feb 14, 2012########################
'''
sripakpa@rescomp-09-171908:/Volumes/Sept_2011_Ext_HD_Parin/minirosetta/Feb_09_SWA_UUAAGU_rebuild_FULL_LENGTH_AND_FINAL_BULGE/main_folder/UUAAGU_Hexaloop$ check_mapper_files_exist.py -N_JOBS 10000 -pattern "FINAL_REBUILD_BULGE/STANDARD/S_&(Process)/mapper_rebuild_bulge.out"
pattern=FINAL_REBUILD_BULGE/STANDARD/S_&(Process)/mapper_rebuild_bulge.out
mapper_outfile (FINAL_REBUILD_BULGE/STANDARD/S_9204/mapper_rebuild_bulge.out) doesn't exist!
mapper_outfile (FINAL_REBUILD_BULGE/STANDARD/S_9538/mapper_rebuild_bulge.out) doesn't exist!
sripakpa@rescomp-09-171908:/Volumes/Sept_2011_Ext_HD_Parin/minirosetta/Feb_09_SWA_UUAAGU_rebuild_FULL_LENGTH_AND_FINAL_BULGE/main_folder/UUAAGU_Hexaloop$ 
'''

##############################9:00 PM Feb 11, 2012###########################
'''
sripakpa@rescomp-09-171908:/Volumes/Sept_2011_Ext_HD_Parin/minirosetta/Feb_09_TRAIL_2_SWA_UUAAGU_rebuild_FULL_LENGTH_AND_FINAL_BULGE/main_folder/UUAAGU_Hexaloop$ check_mapper_files_exist.py -N_JOBS 10000 -pattern "FINAL_REBUILD_BULGE/STANDARD/S_&(Process)/mapper_rebuild_bulge.out" 
pattern=FINAL_REBUILD_BULGE/STANDARD/S_&(Process)/mapper_rebuild_bulge.out
mapper_outfile (FINAL_REBUILD_BULGE/STANDARD/S_9204/mapper_rebuild_bulge.out) doesn't exist!
mapper_outfile (FINAL_REBUILD_BULGE/STANDARD/S_9538/mapper_rebuild_bulge.out) doesn't exist!
sripakpa@rescomp-09-171908:/Volumes/Sept_2011_Ext_HD_Parin/minirosetta/Feb_09_TRAIL_2_SWA_UUAAGU_rebuild_FULL_LENGTH_AND_FINAL_BULGE/main_folder/UUAAGU_Hexaloop$ 
'''
#OK CONSISTENT WITH 3 JOBS RUNNING (INCLUDING REDUCER) BEFORE MANUAL KILL.

#############################################################################

#check_mapper_files_exist.py -N_JOBS 10000 -pattern "FINAL_REBUILD_BULGE/STANDARD/S_&(Process)/mapper_rebuild_bulge.out" > LOG_check_mapper_files.txt


# ~/SWA_RNA_python/SWA_dagman_python/misc/check_mapper_files_exist.py -N_JOBS 1000 -pattern  "REGION_5_4/START_FROM_REGION_5_3/S_&(Process)/region_5_4_sample.out" 

# ~/SWA_RNA_python/SWA_dagman_python/misc/check_mapper_files_exist.py -N_JOBS 1000 -pattern  "DONE/CONDOR/SAMPLER/REGION_5_4_START_FROM_REGION_5_3/&(Process).txt" 


pattern= parse_options( argv, "pattern", "" ) 

N_JOBS= parse_options( argv, "N_JOBS", 0 )

if( N_JOBS<0 ): error_exit_with_message("N_JOBS<0!") 

for queue_ID in range(N_JOBS):

	mapper_outfile=pattern.replace('&(Process)','%d' %(queue_ID) )

	if(exists(mapper_outfile)==False): 
		print "mapper_outfile (%s) doesn't exist!" %(mapper_outfile)
