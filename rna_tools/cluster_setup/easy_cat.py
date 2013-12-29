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

        # make sure to order 0,1,2,... 10, 11, 12, ... 100, 101, ...
        globfiles_with_length = map( lambda x: [len(x),x], globfiles )
        globfiles_with_length.sort()
        globfiles = map( lambda x:x[-1], globfiles_with_length )

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

    cat_file_contributors = cat_file + '.txt'
    fid_txt = open( cat_file_contributors, 'w' )
    for m in range( len( which_files_to_cat[tag] ) ):
        fid_txt.write( '%03d %s\n' % (m,which_files_to_cat[tag][m]) )
    fid_txt.close()

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

