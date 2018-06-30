#!/usr/bin/env python

from sys import argv,exit,stdout,stderr
from os import popen,system
from os.path import exists
import string
import glob

def Help():
    print()
    print()
    print('Usage: cat_outfiles.py <outfile1> <outfile2> ... > <concatenated outfile>')
    print()
    print()
    exit()

if len(argv)<2:
    Help()

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

sequence_line_found    = 0
description_line_found = 0
remark_line_found      = 0
n_file = -1

for out_f in outfiles:
    all_files = glob.glob(out_f)
    for filename in all_files:
        data = open(filename)
        n_file += 1
        for line in data:
            line = line[:-1]
            if not line: break

            if line[:9] == 'SEQUENCE:':
                if sequence_line_found: continue # Should not be any more sequence lines!
                else: sequence_line_found = 1

            if line.find( 'description' ) > -1:
                if description_line_found: continue
                else: description_line_found = 1

            #if line.find( 'REMARK' ) > -1:
            #    if remark_line_found: continue
            #    else: remark_line_found = 1

            description_index = line.find(' S_')
            if description_index < 0:
                description_index = line.find(' F_')

            if description_index >= 0:
                description_index -= 1 # to get rid of space.
                tag = line[description_index:]
                newtag = tag + "_%03d" % n_file
                line = line[:description_index] + newtag

            if len(line) < 1: continue

            fid.write( line+'\n' )


        data.close()

if not final_outfile == "": fid.close()
