#!/usr/bin/python
## look at fa scores

import string
from sys import argv
from os import popen
import sys

file1 = argv[1]

inputfile=open(file1,"r")

firstline=inputfile.readline()
if firstline.find('desc') < 0:
    firstline=inputfile.readline()

if firstline.find('desc') >= 0:
    #Typical out file or score file
    labels=string.split(firstline)
    i=0
    while i < len(labels) :
        print i+1, labels[i]
        i = i +1
else:
   lines = popen( 'grep SCORE '+file1+'| head -n 50').readlines()
   cols = string.split( lines[-1] )
   for i in range( len( cols )/2 ):
       print (2*i+2), cols[2*i]

