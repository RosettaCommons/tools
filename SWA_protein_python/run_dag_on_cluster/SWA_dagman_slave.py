#!/usr/bin/env python

from os import system
from os.path import exists
from sys import argv
from time import sleep
from SWA_util import *

#SWA_dagman_slave.py SLAVE/1/
#SWA_dagman_slave.py SLAVE/2/

jobdir = argv[ 1 ]

assert( exists( jobdir ) )

finished_file = jobdir + '/finished.txt'

# need to do something to show that I'm actually alive!
fid = open( jobdir+'/slave_ok.txt', 'w' )
fid.write( 'OK\n' )
fid.close()

while not exists( finished_file ):

    command_file_name = jobdir + '/run_this_script.txt'

    ran_a_job = 0
    if exists( command_file_name ):

        sleep( 5 ) # stupid file locking.....POTENTIAL BUG HERE..assume that after 5 seconds, the file will be written....

        lines = open( command_file_name ).readlines()
        command = 'source '+lines[0]
        print command
        submit_subprocess( command )
#		if retcode is not 0:
#			print "Rosetta code (c++) terminated by signal" , retcode

        # After done, remove the script file...
        # A more robust signal might be to create a "done" file.
        submit_subprocess( 'rm -rf ' + command_file_name )

#        retcode = system( 'rm -rf ' + command_file_name )
#		if retcode is not 0:
#			print "problem with rm -rf " + command_file_name


        ran_a_job = 1

    sleep( 2 )

