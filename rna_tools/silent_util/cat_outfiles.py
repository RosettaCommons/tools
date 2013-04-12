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

#data = open( outfiles[ 0 ] ).readlines()
#count = 0
#for line in data:
#    fid.write( line )

sequence_line_found    = 0
description_line_found = 0
remark_line_found      = 0

for i in range( len(outfiles) ):

    if not exists( outfiles[i] ):
        stderr.write( 'Does not exist! '+outfiles[i] )
        continue

    data = open(outfiles[i],'r')

    #line = data.readline() # Skip first two lines
    #line = data.readline()
    line = 1

    while line:

        line = data.readline()[:-1]
        if not line: break

        if line[:9] == 'SEQUENCE:':
            if sequence_line_found: continue # Should not be any more sequence lines!
            else: sequence_line_found = 1

        if line.find( 'description' ) > -1:
            if description_line_found: continue
            else: description_line_found = 1

        if line.find( 'REMARK' ) > -1:
            if remark_line_found: continue
            else: remark_line_found = 1

        description_index = line.find('S_')
        if description_index < 0:
            description_index = line.find('F_')
        #if description_index < 0:
        #    description_index = line.find('_')

        if description_index >= 0 and i > 0:
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
