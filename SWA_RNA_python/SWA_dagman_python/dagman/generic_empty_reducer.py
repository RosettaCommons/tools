#!/usr/bin/env python

######################################################################
from SWA_dagman_python.parser.SWA_parse_options import parse_options, parse_seq_num_list_option
from SWA_dagman_python.utility.SWA_util import *

######################################################################
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
######################################################################

start_argv=copy.deepcopy(argv)

done_signal_file = parse_options( argv, "done_signal_file", "" )

if(done_signal_file==""): error_exit_with_message("User need to pass in done_signal_file!")

if(exists(done_signal_file)): error_exit_with_message("done_signal_file (%s) already exist!" %(done_signal_file))


create_generic_done_signal_file(done_signal_file)


print "----------------------------------------------------------------------------------------------------------------------------"
print "Successfully RAN %s" %( list_to_string(start_argv) )
print "----------------------------------------------------------------------------------------------------------------------------"

