#!/usr/bin/env phenix.python
import os.path
import imp

file_path = os.path.split( os.path.abspath(__file__) ) [0]
imp.load_source('erraser_wrapper', file_path + '/erraser_wrapper.py')
imp.load_source('erraser_option', file_path + '/erraser_option.py')
imp.load_source('erraser_util', file_path + '/erraser_util.py')

from erraser_wrapper import erraser
from erraser_option import erraser_option
from erraser_util import *

option = erraser_option()
option.read_cmdline_erraser( sys.argv )
erraser( option )
