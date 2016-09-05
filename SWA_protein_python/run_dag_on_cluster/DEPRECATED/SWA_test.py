#!/usr/bin/env python

from sys import argv,exit,stdout, stderr
import sys
import string
#from os import system,popen
import subprocess
from os.path import basename,exists
from random import randrange
from time import sleep
from os import system,popen2,getcwd
from os import system,popen
from SWA_util import *
from SWA_filter_silent_file import *
from SWA_cat_outfiles import *
from glob import glob


cat_outfile="REGION_1_6/start_from_region_2_6_sample.out"
filter_outfile="laptop_test.txt"
score_diff_cut=8.0

(actual_num_structures, actual_score_gap)=filter_silent_file(silent_file=cat_outfile, outfile_name=filter_outfile, max_score_gap=score_diff_cut)
print "actual_num_structures=%d, actual_score_gap=%f " %(actual_num_structures, actual_score_gap)


#/home/sripakpa/SWA_dagman_python/SWA_sampling_post_process.py -indir_prefix REGION_2_1_test/START_FROM_REGION_2_0 -score_diff_cut   24.000


#filter_silent_file(silent_file='/home/parin/minirosetta/Mar_7_2r8s_suite_cluster_loop_007_suite_010/REGION_5_2/start_from_region_5_1_sample.out', outfile_name='test.txt', max_score_gap=10)

#(actual_num_structures, actual_score_gap)=filter_silent_file(silent_file='REGION_5_2/start_from_region_5_1_sample.out', outfile_name='REGION_5_2/test_1.txt', max_n_struct=40000)

#print "actual_num_structures=%d, actual_score_gap=%f " %(actual_num_structures, actual_score_gap)

#num_struct=filter_silent_file(silent_file='REGION_2_0/START_FROM_REGION_0_0/region_2_0_sample.out', outfile_name='REGION_2_0/test_1.txt', max_n_struct=40000)
#print "num_struct=", num_struct


#concatenate_outfiles(['','region_0_2_sample.cluster.out'], 'test.txt')

'''
silent_file= 
outfile_name=
REVERSE=  True
scorecol_name=score
mode_max_struct, max_n_struct=40000, max_score_gap=0
score is located at column_num: 2 tag is located at column_num: 23
'''

#SCORE_silentfile="/home/parin/minirosetta/Mar_7_2r8s_suite_cluster_loop_007_suite_010/REGION_5_2/SCORE_start_from_region_5_1_sample.out"

#SCORE_silentfile="SCORE_start_from_region_5_1_sample.out"


#open( SCORE_silentfile, 'r' )

#filter_silent_file(silent_file, outfile_name, scorecol_name="score", REVERSE=True, max_n_struct=0, max_score_gap=0):


'''
lines = popen("tail -n 20 0.out").readlines()
print lines		

for line in lines:
#	line=line[:-1]
	if (line == 'Exit job, please check errfile to ensure that there are no errors.' ):  #There must be a better way to write this!
		print "completed!"                  
		break
'''
