#!/usr/bin/python


import sys
import string
from os import system,popen
from os.path import basename, abspath, dirname
from glob import glob

outfiles = sys.argv[1:]

scripts_path = dirname( abspath( sys.argv[0] ) )

which_files_to_cat = {}

for outfile in outfiles:
    if (outfile[-4:] == '.out' ):
        #Old style, user specified a bunch of outfiles.
        tag = string.join( string.split( outfile,'_' )[:-2] , '_')
        if tag not in which_files_to_cat.keys():
            which_files_to_cat[tag] = []
        which_files_to_cat[tag].append( outfile )

    else: #New style -- look inside a directory!

        globfiles = glob( outfile+'/*/*out' )

        #print globfiles
        if len( globfiles ) == 0: globfiles = glob( outfile + '/*out'  )
        #print globfiles

        globfiles.sort()

        for file in globfiles:
            tag = basename( file ).replace('.out','')
            if tag not in which_files_to_cat.keys():
                which_files_to_cat[tag] = []
            which_files_to_cat[tag].append( file )


for tag in which_files_to_cat.keys():
    #if (len( which_files_to_cat[tag] ) == 1) : continue

    cat_file = tag+".out"
    #    cat_file = "cat_"+tag+".out"
    print "Catting into: ",cat_file,
    command = 'python %s/cat_outfiles.py %s >  %s ' % \
              (scripts_path, string.join( which_files_to_cat[tag] ) ,
               cat_file )
    #print command
    system( command )

    lines = popen( 'grep SCORE '+cat_file).readlines()
    print '... from %d primary files. Found %d  decoys.' % (len( which_files_to_cat[tag] ),len(lines)-1)

    fid_sc = open( cat_file.replace('.out','.sc'),'w' )
    for line in lines:
        fid_sc.write( line )
    fid_sc.close()

for outfile in outfiles:
    command = 'rm -rf '+outfile
    #print command
    #system( command )

