#!/usr/bin/env python

from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
from SWA_parse_options import parse_options
from SWA_util import *

print "print len(argv)=", len(argv)

print
print

silent_file = parse_options( argv, "silent_file", "" )
assert( exists( silent_file ) )
remove_variant_types = parse_options( argv, "remove_variant_types", "true" )

cluster_rmsd= parse_options( argv, "cluster_rmsd", 2.0 )

common_args=""
common_args_file = parse_options(argv, "common_args", "")
if(common_args_file != ""):
	print "common_args_file= ", common_args_file
	common_args=open( common_args_file  ).readlines()[0][:-1]
#	common_args=popen_and_readlines(common_args_location)[0][:-1]
	
print common_args

###After parsing there should only be one argv left which is the python file name: SWA_rna_build_dagman.py####
if(len(argv)!=1):
	print argv
	print "leftover len(argv)=", len(argv)
	assert(False)

print
print

#Optional



cluster_rmsd_int=int(cluster_rmsd*100)

output_filename='cluster_%d_%s' % (cluster_rmsd_int , basename(silent_file))

print "%f Angstrom" % cluster_rmsd

system( 'rm %s'  % output_filename )

bin_folder="~/src/mini/bin/"
compile_type="release" #"debug" # 
rna_swa_test= bin_folder + "/rna_swa_test.linuxgcc" + compile_type

command=rna_swa_test
command += " -algorithm cluster_old"
command += " -in:file:silent " + silent_file
command += " -in:file:silent_struct_type  binary_rna"
command += " -database ~/minirosetta_database"
command += " -out:file:silent " + output_filename
command += " -score_diff_cut 1000000.000"

if(common_args==""):
	command += " -cluster:radius %s " % cluster_rmsd
else:
	command += " -loop_cluster_radius %s " % cluster_rmsd
	command += " %s " %common_args

command += " -nstruct 1000"

command += ' > cluster_output.txt'
print command 
submit_subprocess( command )

print 'SWA_extract_pdb.py -silent_file %s -remove_variant_types %s' % (output_filename, remove_variant_types)
submit_subprocess( 'SWA_extract_pdb.py -silent_file %s -remove_variant_types %s' % (output_filename, remove_variant_types) )

