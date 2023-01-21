#!/usr/bin/env python
## look at fa scores

import string
from sys import argv
from os import popen
import sys

def Help():
    print( argv[0] + ' <scorefilename>' )
    print( )
    print( '  Used for Rosetta silent files or score files -- figure out' )
    print( '  which column numbers correspond to which score terms.' )
    print( )
    exit( 0 )

if len(argv) < 2:
    Help()

file1 = argv[1]

inputfile=open(file1,"r")

firstline=inputfile.readline()
if firstline.find('desc') < 0:
    firstline=inputfile.readline()

if firstline.find('desc') >= 0:
    #Typical out file or score file
    labels=firstline.split()
    i=0
    while i < len(labels) :
        print( '%4d %s' % ( i+1, labels[i] ) )
        i = i+1
else:
    lines = popen( 'grep SCORE '+file1+'| head -n 50').readlines()
    if len( lines ) > 0:
        cols =  lines[-1].split()
        for i in range( len( cols )/2 ):
            print( '%4d %s' %  (2*i+2, cols[2*i] ) )
    else:
        lines = popen( 'head -n 1 '+file1 ).readlines()
        cols =  lines[0].split()
        for i in range( len( cols ) ):
            print( '%4d %s' % (i+1, cols[i]) )


print( )
print( 'Score tags for: ',file1 )
