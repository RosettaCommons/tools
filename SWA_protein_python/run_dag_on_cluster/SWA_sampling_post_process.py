#!/usr/bin/env python

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
from os.path import basename, dirname, exists, expanduser
from time import sleep

from SWA_util import *
from SWA_cat_outfiles import *
from SWA_filter_silent_file import *
from SWA_parse_options import parse_options

#############################################################################

python_command=get_python_command(argv)
print_title_text("Enter " + python_command)

score_diff_cut = parse_options( argv, "score_diff_cut", 0.0 )
print "score_diff_cut= %8.3f" %score_diff_cut

if(score_diff_cut<0.1):
	error_exit_with_message("score_diff_cut (%8.3f) < 0.1 !" %(score_diff_cut))


indir_prefix = parse_options( argv, "indir_prefix", "" )
if(indir_prefix == ""):
	error_exit_with_message("need to input indir_prefix!")


if(len(argv)!=1):
	print argv
	print "leftover len(argv)=", len(argv)
	assert(False)

#indir_prefix = argv[1]




########################################################################3
globstring = indir_prefix+'*/*sample.out'
print "globstring= ", globstring
sys.stdout.flush()

globfiles = glob( globstring )
if len( globfiles ) == 0:
	sleep( 10 ) #IN what situation is this needed
	globfiles = glob( globstring )
globfiles.sort()

if len( globfiles ) == 0: ##OK still no outfiles....probably a bug with rosetta output
	error_exit_with_message("len(globfiles)==0 for globstring %s !!" %globstring)

#print globfiles

cat_outfile = dirname(indir_prefix) + '/' + basename(indir_prefix).lower() + '_sample.out'


if exists( cat_outfile ):
	submit_subprocess("rm -r %s" %(cat_outfile))	

PYDIR = expanduser('~')+'/SWA_dagman_python/'              #######################################CRAP THIS IS BUGGY AGAIN>>>>>>###############################
assert( exists( PYDIR ) )


#command = PYDIR+'/SWA_cat_outfiles.py '+string.join( globfiles ) + ' -o ' + cat_outfile  
#print( command )
#print("cat_outfile=%s " %(cat_outfile) )

concatenate_outfiles(infile_list=globfiles, outfile=cat_outfile)

##########################################
##########################################
# We don't need the whole file -- just
filter_outfile = cat_outfile.replace('.out','_filtered.out')  ####BEWARE THIS IS SPECIFIC FORM IS HARD CODED in SWA_rna_build_dagman.py

#wait_for_file( cat_outfile )
#command = PYDIR+'/SWA_filter_silent_file.py'
#command +=' -silent_file %s' %(cat_outfile)
#command +=' -outfile_name %s' %(filter_outfile)
#command +=' -n_struct 40000'   ####Allow score_gap mode in the future
#print command
#submit_subprocess(command)

(actual_num_structures, actual_score_gap)=filter_silent_file(silent_file=cat_outfile, outfile_name=filter_outfile, max_score_gap=score_diff_cut)
print "actual_num_structures=%d, actual_score_gap=%f " %(actual_num_structures, actual_score_gap)

if(actual_num_structures < 8001):
	print "Too few structure kept!, refilter_silent_file"
	(actual_num_structures, actual_score_gap)=filter_silent_file(silent_file=cat_outfile, outfile_name=filter_outfile, max_n_struct=8000)
	print "Refiltered: actual_num_structures=%d, actual_score_gap=%f " %(actual_num_structures, actual_score_gap)

sys.stdout.flush()
##########################################

# DISK SPACE!

###Mod this for testing for now

####Delete all the output folders#######

globstring = indir_prefix+'*'
globfiles = glob( globstring )
globfiles.sort()
for globfile in globfiles:
    command = 'rm -rf '+globfile
    print( command )
    submit_subprocess( command ) 


###Delete the cat_outfile#############

command = 'rm '+cat_outfile
print( command )
submit_subprocess( command )


print_title_text("Exit " + python_command)

