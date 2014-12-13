#!/usr/bin/env python

######################################################################
from SWA_dagman_python.utility.SWA_util import *


######################################################################

print "Enter simple_test_script.py"

sys.stdout.flush()
sys.stderr.flush()

system('echo "BLAH" > simple_test_script_echo.txt') 

print "CURRENT_DIRECTORY=%s" %(os.path.abspath("."))

sys.stdout.flush()
sys.stderr.flush()

sleep_time=60

print "sleep %s seconds" %(sleep_time)

sys.stdout.flush()
sys.stderr.flush()

sleep(sleep_time)

print "Exit simple_test_script.py"

sys.stdout.flush()
sys.stderr.flush()

