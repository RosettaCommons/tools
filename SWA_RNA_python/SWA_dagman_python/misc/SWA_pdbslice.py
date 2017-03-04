#!/usr/bin/env python

from os import popen,system
from sys import argv,stderr
import string
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################

############################
if(argv.count('-segments')==0): error_exit_with_message("argv.count('-segments')==0")

pos = argv.index('-segments')
del argv[pos]

subset_residues = []

segment_residues = []

goodint = 1
while goodint:
    try:
        segment_residue = int(argv[pos])
        segment_residues.append( segment_residue )
        del argv[pos]
        #stderr.write('%d ' % segment_residue )
    except:
        goodint = 0

for i in range( len(segment_residues)/2):
	for j in range( segment_residues[2*i], segment_residues[2*i+1]+1 ):
		subset_residues.append( j )

print 'slice_segments: ', subset_residues

prefix= parse_options( argv, "prefix", "" ) 

if(prefix==""):

	pdbfiles = argv[1:-1]

	if(len(pdbfiles)!=1): error_exit_with_message("prefix=="" but len(pdbfiles)=(%s)!=1" %(len(pdbfiles)))

	output_pdb = argv[-1]

else:

	pdbfiles = argv[1:]



for pdbfile in pdbfiles:

	######################################
	if(prefix):
		outid = open(prefix+pdbfile,'w')
	else:  
		outid = open(output_pdb,'w')
	#######################################

	if(pdbfile[-4:] != '.pdb'): error_exit_with_message("pdbfile[-4:] != '.pdb'")

	lines = open(pdbfile).readlines()

	i = 0
	old_seqnum = '   '
	for line in lines:
		if(len(line)<30): error_exit_with_message("len(line)<30")

		curr_seqnum = line[22:26]

		if(curr_seqnum != old_seqnum):
			i += 1
			old_seqnum = curr_seqnum

		if( i not in subset_residues ): continue

		try:
			curr_seqnum_int = int(curr_seqnum)
		except:
			error_exit_with_message("Cannot convert curr_seqnum (%s) to int!" %(curr_seqnum))

		if( (curr_seqnum_int < 1) or (curr_seqnum_int > 10000000000) ): 
			error_exit_with_message("(curr_seqnum_int) < 1 or (curr_seqnum_int > 10000000000) | curr_seqnum=%s" %(curr_seqnum))

		outid.write(line)

	outid.close()

