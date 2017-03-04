#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
import os
import copy
from glob import glob
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################

#SWA_pdb_slice_wrapper.py -pdb_pattern DAS_FARFAR_FULL_^^.pdb -segments 1 13 20 29


pdb_pattern = parse_options( argv, "pdb_pattern", "" )

segments=parse_options( argv, "segments", [0] )

pdb_pattern=pdb_pattern.replace("^^","*")


submit_subprocess("copy_AL_Hahimi_models.py")

print "segments= %s " %(segments)
print "pdb_pattern= %s " %(pdb_pattern)

if( len(segments)%2 != 0): error_exit_with_message(len(segments)%2 != 0)

pdb_list = glob( pdb_pattern )	
		
pdb_list.sort()



for pdb in pdb_list:

	output_pdb=pdb.replace("FULL", "TAR_BULGE")

	pdbslice_command="SWA_pdbslice.py %s -segments %s %s" %(pdb, list_to_string(segments), output_pdb) 
	submit_subprocess(pdbslice_command)

	submit_subprocess("renumber_pdb_in_place.py %s" %(output_pdb))
