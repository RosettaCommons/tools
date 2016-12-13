#!/usr/bin/env python

######################################################################
from SWA_DAG_util import *

######################################################################

from SWA_dagman_python.misc.SWA_cat_outfiles import *
from SWA_dagman_python.misc.SWA_filter_silent_file import *
from SWA_dagman_python.utility.SWA_util import *
######################################################################

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
from os.path import basename, dirname, exists, expanduser
from time import sleep


######################################################################
######Feb 12, 2012: Should switch from using indir_prefix to obtaining the equivalent information from the condor_submit_file!
def delete_sampler_outfiles_and_folders(indir_prefix, delete_files):

	#Updated on Jan 10, 2012 (change '*' to '_S_*/'), this prevent possible confusion of START_FROM_REGION_0_0_AND_2_0_S_100/ with START_FROM_REGION_0_0_S_100/
	#Also to avoid possibility of deleting the cat_outfile and filter_outfile which have similar name to outfolder.

	'''
	if(delete_files):

		folder_globstring = indir_prefix+'/'
		folder_globfiles = glob( folder_globstring )
		folder_globfiles.sort()
		for folder_globfile in folder_globfiles:
				command = 'rm -rf '+ folder_globfile
				print( command )
				sys.stdout.flush()
				sys.stderr.flush()
				submit_subprocess( command )
	'''

	if(delete_files):
		outfolder=indir_prefix +'/'
		#if(exists(outfolder)==False): error_exit_with_message("outfolder (%s) doesn't exist!" %(outfolder))
		if( not exists(outfolder) ):
			print "WARNING, outfolder (%s) doesn't exist!" % (outfolder)
			return False
		print "rm -r %s" %(outfolder)
		submit_subprocess("rm -r %s" %(outfolder))
		return True

#############################################################################
def concatenate_sampler_silent_file(condor_submit_file, cat_outfile):

	#######REWROTE THIS FUNCTION ON FEB 08, 2012##########

	if(exists( cat_outfile )): error_exit_with_message("cat_outfile (%s) already exist!" %(cat_outfile) )

	########################################################
	#less COMBINE_DS_REGIONS//COMBINE_REGION_8_0_and_1_7//BY_APPEND/S_1/combine_region_8_0_and_1_7_sample.out
	#Error!: silent_file (COMBINE_DS_REGIONS//COMBINE_REGION_8_0_and_1_7//BY_APPEND/S_0/combine_ds_regions_sample.out) doesn't exist!

	if(exists(condor_submit_file)==False): error_exit_with_message("condor_submit_file (%s) doesn't exist!" %(condor_submit_file) )

	lines = safe_readlines(condor_submit_file)

	num_queue_line=0
	num_mapper_outfiles_line=0

	num_mapper_jobs=0
	mapper_outfile=""

	for line in lines:

		cols=line.split()

		if(len(cols)<2): continue # error_exit_with_message("len(cols)<2 for line (%s)" %(line))

		if(cols[0] == 'Queue'):

			num_queue_line+=1

			if( len(cols)!=2 ): error_exit_with_message("len(cols)!=2 for Queue line=(%s)" %(line))

			num_mapper_jobs=int(cols[1])

		if( check_tag( cols[0], "mapper_outfiles")  ):

			num_mapper_outfiles_line=+1

			#For general case, can have more than one mapper_outfile, but for SAMPLER and DS_COMBINE_SAMPLER, only have 1 mapper_outfile [per mapper_job].
			if( len(cols)!=3 ): error_exit_with_message("len(cols)!=3 for mapper_outfiles line=(%s)" %(line))

			if( cols[1] != "=" ): kill_all_slave_jobs_and_exit("cols[1] != \"=\" for mapper_outfiles line (%s)" %(line))

			mapper_outfile  = cols[2]


	if(num_queue_line!=1): error_exit_with_message("num_queue_line=(%s)!=1" %(num_queue_line))
	if(num_mapper_outfiles_line!=1): error_exit_with_message("num_mapper_outfiles_line=(%s)!=1" %(num_mapper_outfiles_line))

	if(mapper_outfile==""): error_exit_with_message("mapper_outfile=\"\"")
	if(num_queue_line<=0): error_exit_with_message("num_queue_line<=0")

	####SOME POSSIBLE mapper_outfiles examples######
	#mapper_outfiles =  REGION_0_1/START_FROM_REGION_0_0/ONLY_JOB/region_0_1_sample.out

	#mapper_outfiles =  REGION_1_4/START_FROM_REGION_2_4/S_$(Process)/region_1_4_sample.out

	#mapper_outfiles =  COMBINE_DS_REGIONS//COMBINE_REGION_6_2_and_3_5//BY_PREPEND/S_$(Process)/combine_region_6_2_and_3_5_sample.out


	########################################################
	non_empty_silent_file_list=[]

	for queue_ID in range(num_mapper_jobs):

		silent_file=mapper_outfile.replace('$(Process)','%d' %(queue_ID) )

		if(exists(silent_file)==False): error_exit_with_message("silent_file (%s) doesn't exist!" %(silent_file))

		silent_data=safe_open(silent_file, mode='r', Is_master=False)

		first_silent_line=silent_data.readline()

		silent_data.close()

		Is_empty=False

		#if first_silent_line.index( "empty silent_file" ) > -1:
		if "empty silent_file" in first_silent_line:
			Is_empty=True
		else:
			assert_is_valid_non_empty_silent_file(silent_file)
			non_empty_silent_file_list.append(silent_file)

		print "queue_ID=%d | silent_file=%s | Is_empty=%s" %(queue_ID , silent_file, Is_empty)

	print "num_mapper_jobs=%d | len(non_empty_silent_file_list)=%d" %(num_mapper_jobs, len(non_empty_silent_file_list))

	sys.stdout.flush()
	sys.stderr.flush()

	if(len(non_empty_silent_file_list)==0): return True #Indicate that there are no non_empty silent_file

	concatenate_outfiles(infile_list=non_empty_silent_file_list, outfile=cat_outfile)

	return False #Indicate that there are non-empty silent_files!


	'''##PRE FEB 08, 2012 VERSION
	########################################################################3
	globstring = indir_prefix+'/*/*sample.out'  #Updated on Jan 10, 2012
	print "globstring= ", globstring
	sys.stdout.flush()

	globfiles = glob( globstring )
	if len( globfiles ) == 0:
		sleep( 10 ) #IN what situation is this needed
		globfiles = glob( globstring )
	globfiles.sort()

	####################################################################################################################################
	if(len( globfiles ) == 0): ##OK still no outfiles....probably a bug with rosetta output
		return True

	####################################################################################################################################
	for i in range(0, len(globfiles)): #Output the files!

		print "i=%d, globfiles[i]=%s" %(i , globfiles[i])

		######Sept 24, 2011 CONSISTENCY CHECK##########
		if(globfiles[i].count("DEBUG")!=0): error_exit_with_message( 'globfiles[%d].count("DEBUG")!=0, globfiles[%d]=%s' %(i, i, globfiles[i]) )
	####################################################################################################################################

	sys.stdout.flush()
	sys.stderr.flush()

	if(exists( cat_outfile )):
		error_exit_with_message("cat_outfile (%s) already exist!" %(cat_outfile) )
		#print "Warning cat_outfile (%s) already exist......removing!" %(cat_outfile)
		#submit_subprocess("rm -r %s" %(cat_outfile))

	concatenate_outfiles(infile_list=globfiles, outfile=cat_outfile)

	return False #zero_sampler_file
	'''

#############################################################################

def filter_sampler_cat_outfile(filter_outfile, cat_outfile, min_filtered_struct):

	remove_SCORE_file=False #Might want to keep this False while testing!

	if(exists(filter_outfile)):
		print "Warning filter_outfile (%s) already exist....removing!" %(filter_outfile)
		submit_subprocess("rm %s" %(filter_outfile))

	#Remove all score_filtering on Jan 12, 2012.

	(actual_num_structures, actual_score_gap)=filter_silent_file(silent_file=cat_outfile, outfile_name=filter_outfile, max_n_struct=min_filtered_struct, remove_SCORE_file=remove_SCORE_file)
	print "In filtered_silent_file: actual_num_structures=%d, actual_score_gap=%f " %(actual_num_structures, actual_score_gap)

	sys.stdout.flush()
