#!/usr/bin/python

from glob import glob
from os import getcwd,chdir,system
from os.path import basename,dirname

outdir_name = 'SWM'
globdirs = glob( '*/%s' % outdir_name )
CWD = getcwd()
for globdir in globdirs:
    chdir( dirname( globdir ) )
    print dirname( globdir )
    system( 'easy_cat.py %s' % outdir_name )
    print
    chdir( CWD )
