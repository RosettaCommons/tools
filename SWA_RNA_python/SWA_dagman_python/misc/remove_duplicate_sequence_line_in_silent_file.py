#!/usr/bin/env python

#from scipy import savetxt,loadtxt, array
from sys import argv,exit
from os import popen, system
from os.path import basename,exists,expanduser
import string
import commands
from glob import glob
######################################################################

assert( len(argv)>1)
file_name = argv[1]
file_name_correct=file_name + '.correct'
system( 'rm %s' % file_name_correct)


#head_file = file_name + '.head10'
#system( 'head -n 10 %s > %s' % (file_name, head_file) )

infile = open(file_name, 'r')

seq_line=infile.readline()

infile.close()

print seq_line, len (seq_line)


outfile = open(file_name_correct, 'w')
outfile.write(seq_line);
  
#print 'grep -v %s %s > %s' % (seq_line, file_name, file_name_correct)

system( 'grep -v "%s" %s >> %s' % (seq_line, file_name, file_name_correct))

system( 'rm -r %s' % file_name)
system( 'mv %s %s' % (file_name_correct, file_name))
