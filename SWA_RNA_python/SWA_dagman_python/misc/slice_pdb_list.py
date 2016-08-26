#!/usr/bin/env python

from os import system,popen
from sys import argv,exit
import sys
import traceback
from glob import glob
import string
from os.path import basename, dirname, exists, expanduser
from time import sleep


globstring ='C*.pdb'
print "globstring= ", globstring


globfiles = glob( globstring )
if len( globfiles ) == 0:
	sleep( 10 ) #IN what situation is this needed
	globfiles = glob( globstring )
globfiles.sort()



print globfiles

for pdb in globfiles:
	command= 'pdbslice.py %s -segments 3 7 12 16 trimmed_' %(pdb)
	system(command)
	command= 'renumber_pdb_in_place.py  trimed_%s' %(pdb)

