#!/usr/bin/env python

import string
from sys import argv
from os.path import exists,basename

from SWA_util import *

###Example########
###/home/sripakpa/SWA_dagman_python/SWA_sampling_pre_process.py REGION_5_1 REGION_5_0 CONDOR/REGION_5_1_START_FROM_REGION_5_0.condor

####WARNING#######
#This code is extremely not robust..Will rewrite this when I have the time. Mar 21, 2010

python_command = get_python_command(argv)
print_title_text("Enter: " + python_command)


outdir = argv[1]  
dir_prev = argv[2] 
condor_submit_file = argv[3] 

# Need directories to be setup. Could do this with a preprocessing script? That would also allow
# for convenient setup of Queue number in the condor script.

#lines = popen( 'grep SCORE %s_sample.cluster.out' % ( dir_prev.lower() ) ).readlines()
#lines = popen_and_readlines( 'grep SCORE %s_sample.cluster.out' % ( dir_prev.lower() ) )

#################################################################################################################
input_silent_file='%s_sample.cluster.out' %(dir_prev.lower())  ###AGAIN this is hard coded 

firstlines = popen_and_readlines('head -n 3 '+ input_silent_file)

col_name_line=firstlines[1] ###AGAIN this assume that the second line of the silent file is the column name line
col_name_list=string.split(col_name_line)

try:
	tag_col_index=col_name_list.index('description')
except:
	error_exit_with_message("Cannot find description column index! ")

print "tag is located at column_num: %d" %(tag_col_index+1)

SCORE_silentfile='SCORE_' + basename(input_silent_file)

command = 'grep SCORE %s > %s' %(input_silent_file, SCORE_silentfile)

submit_subprocess(command)

try:
	infile = open( SCORE_silentfile,'r')
except:
	error_exit_with_message("cannot open %s " % SCORE_silentfile)
	
###Note: Elsewhere, the python code had been hard coded to assume that the tags from input_silentfile is of the form:
### S_0, S_1, S_2, S_3,...... Personally I think we should relax this assumption but  to do thiswill have to make 
### many changes elsewhere in the code. The best I could do for now is to ensure here that the input_silentfile indeed
### does have this specific form and raise an error otherwise.

pose_num=-1

for line in infile:
	
	if(Is_column_name_line(line)==True): continue 

	pose_num+=1
	assumed_tag = 'S_%d' %pose_num
#	assumed_tag = 'S_%06d' %pose_num

#S_000000

	cols = string.split( line )
#	print cols

	actual_tag=cols[tag_col_index]

	if(assumed_tag==actual_tag):
		newdir = outdir+'/START_FROM_%s_%s' % ( dir_prev.upper(), actual_tag )   #ONE NAME CHANGE and the code will break completely...
		if not exists( newdir ):  submit_subprocess( 'mkdir -p '+newdir )
	else:
		error_exit_with_message("assumed_tag (%s)!=actual_tag (%s)" %(assumed_tag, actual_tag) )	

infile.close()
submit_subprocess("rm %s" %SCORE_silentfile)


##################################################################################################################
N_JOBS = pose_num+1

print "N_JOBS (# pose from previous job cluster silent_file): %d" %N_JOBS

# Go through condor submission file and update number of jobs...
# May need to be careful about any special cases.
lines = open( condor_submit_file ).readlines()
fid = open( condor_submit_file, 'w' ) #delete the old file and create a new blank one to be written.

for line in lines:
    if len( line ) > 5 and line[:5] == 'Queue':
        fid.write( 'Queue %d\n' % N_JOBS )
    else:
        fid.write( line )

fid.close()

print_title_text("Exit: " + python_command)
