#!/usr/bin/env python
import sys
import string
from os import system,popen,getcwd,chdir,listdir
from os.path import basename,abspath, dirname,exists,join,isdir
from glob import glob

outfiles = sys.argv[1:]
scripts_path = dirname( abspath( sys.argv[0] ) )

which_files_to_cat = {}

for outfile in outfiles:
    if not exists( outfile ): # look inside subdirectories
        d = '.'
        CWD = getcwd()
        subdirs = [o for o in listdir(d) if isdir(join(d,o))]
        for subdir in subdirs:
            chdir( subdir )
            print(subdir)
            system( 'easy_cat.py %s' % outfile )
            print()
            chdir( CWD )
        continue
    if (outfile[-4:] == '.out' ):
        #Old style, user specified a bunch of outfiles.
        tag = "_".join(outfile.split('_')[:-2])
        if tag not in which_files_to_cat.keys():
            which_files_to_cat[tag] = []
        which_files_to_cat[tag].append( outfile )

    else: #New style -- look inside a directory!

        globfiles = glob( outfile+'/*/*out' )

        #print globfiles
        if len( globfiles ) == 0: globfiles = glob( outfile + '/*out'  )
        #print globfiles

        # Remove any "checkpoint" files from stepwise checkpointing
        globfiles = filter( lambda x: "S_" not in x and "_checkpoint" not in x, globfiles )

        # make sure to order 0,1,2,... 10, 11, 12, ... 100, 101, ...
        globfiles_with_length = list(map( lambda x: [len(x),x], globfiles ))
        globfiles_with_length.sort()
        globfiles = list(map( lambda x:x[-1], globfiles_with_length ))

        for file in globfiles:
            tag = basename( file ).replace('.out','')
            if tag not in which_files_to_cat.keys():
                which_files_to_cat[tag] = []
            which_files_to_cat[tag].append( file )


for tag in which_files_to_cat.keys():
    #if (len( which_files_to_cat[tag] ) == 1) : continue

    cat_file = tag+".out"
    print("Catting into: ",cat_file, end='')
    command = 'python %s/cat_outfiles.py %s >  %s ' % \
              (scripts_path, " ".join( which_files_to_cat[tag] ) ,
               cat_file )
    system( command )

    cat_file_contributors = cat_file + '.txt'
    fid_txt = open( cat_file_contributors, 'w' )
    for m in range( len( which_files_to_cat[tag] ) ):
        fid_txt.write( '%03d %s\n' % (m,which_files_to_cat[tag][m]) )
    fid_txt.close()

    lines = popen( 'grep SCORE '+cat_file).readlines()
    print('... from %d primary files. Found %d  decoys.' % (len( which_files_to_cat[tag] ),len(lines)-1))

    fid_sc = open( cat_file.replace('.out','.sc'),'w' )
    for line in lines:
        fid_sc.write( line )
    fid_sc.close()

for outfile in outfiles:
    command = 'rm -rf '+outfile
    #print command
    #system( command )

