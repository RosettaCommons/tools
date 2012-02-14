#!/usr/bin/python

from sys import argv,exit,stdout,stderr
from os import popen,system
from os.path import exists
import string

def Help():
    print
    print
    print 'Usage: cat_outfiles.py <outfile1> <outfile2> ... > <concatenated outfile>'
    print
    print
    exit()

#if len(argv)<3:
#    Help()

outfiles = argv[1:]

final_outfile = ""
if outfiles.count( "-o" ) > 0:
    pos = outfiles.index( "-o" )
    del( outfiles[ pos ] )
    final_outfile = outfiles[ pos ]
    del( outfiles[ pos ] )

if final_outfile == "":
    fid = stdout
else:
    fid = open( final_outfile, "w" )

data = open( outfiles[ 0 ] ).readlines()
for line in data:
    fid.write( line )


for i in range(1, len(outfiles)):

    if not exists( outfiles[i] ):
        stderr.write( 'Does not exist! '+outfiles[i] )
        continue

    data = open(outfiles[i],'r')

    line = data.readline() # Skip first two lines
    line = data.readline()

    while line:
        line = data.readline()[:-1]

        if line[:9] == 'SEQUENCE:': continue # Should not be any more sequence lines!

        description_index = line.find('S_')
        if description_index < 0:
            description_index = line.find('F_')
        #if description_index < 0:
        #    description_index = line.find('_')


        if description_index >= 0:
            tag = line[description_index:]

            tagcols = string.split(tag,'_')
            try:
                tagnum = int( tagcols[-1] )
                tagcols[-1] = '%04d_%06d' %  ( i, tagnum )
                newtag = string.join( tagcols,'_')

                line = line[:description_index] + newtag
            except:
                continue

        if len(line) < 1: continue

        fid.write( line+'\n' )

    data.close()

if not final_outfile == "": fid.close()
