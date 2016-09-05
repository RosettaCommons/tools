#!/usr/bin/env python

from sys import argv,exit
import sys
import traceback
import string
import os
from os.path import basename, dirname, exists, expanduser
from time import sleep
import copy
import string
from os import popen 
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
from SWA_dagman_python.utility.Rosetta_to_standard_PDB_functions import *
############################################################################################


input_pdb_file= parse_options( argv, "s", "" )

remove_hydrogen= parse_options( argv, "remove_hydrogen", "False" )

output_pdb_file= parse_options( argv, "output_pdb", "")

Rosetta_to_standard_PDB_func(input_pdb_file, remove_hydrogen, output_pdb_file, VERBOSE=True)
############################################################################################





########ATOM     41  O4*  rG A   2     -10.469  -6.154   2.097  1.00 10.00              
########ATOM     24 1H5*  rG A   1     -12.562  -6.584   9.293  1.00  0.00  
########12345678901234567890123456789012345678901234567890123456789012345678901234567890
########0        1         2         3         4         5         6         7         8





'''
01234567890123456789012345678901234567890123456789012345678901234567890123456789	
ATOM      1  P    rG R   1     101.133 -43.649  16.700  1.00 59.94           P  
ATOM      2  O1P  rG R   1     101.879 -44.726  16.006  1.00 60.94           O  
ATOM      3  O2P  rG R   1     101.834 -42.388  17.023  1.00 56.97           O  
ATOM      4  O5*  rG R   1      99.855 -43.236  15.828  1.00 56.83           O  
ATOM      5  C5*  rG R   1      99.271 -44.139  14.899  1.00 55.00           C  
ATOM      6  C4*  rG R   1      98.207 -43.454  14.054  1.00 53.51           C  
ATOM      7  O4*  rG R   1      97.351 -42.626  14.888  1.00 51.34           O  
ATOM      8  C3*  rG R   1      98.731 -42.485  13.000  1.00 52.08           C  
ATOM      9  O3*  rG R   1      99.328 -43.159  11.846  1.00 51.13           O  
ATOM     10  C2*  rG R   1      97.455 -41.706  12.689  1.00 51.48           C  
ATOM     11  O2*  rG R   1      96.588 -42.362  11.784  1.00 51.16           O  
ATOM     12  C1*  rG R   1      96.835 -41.578  14.089  1.00 52.18           C  
ATOM     13  N9   rG R   1      97.175 -40.288  14.674  1.00 52.27           N  
ATOM     14  C8   rG R   1      98.138 -40.006  15.612  1.00 52.60           C  
ATOM     15  N7   rG R   1      98.207 -38.733  15.905  1.00 54.08           N  
ATOM     16  C5   rG R   1      97.239 -38.137  15.109  1.00 51.50           C  
ATOM     17  C6   rG R   1      96.861 -36.780  14.975  1.00 51.23           C  
ATOM     18  O6   rG R   1      97.309 -35.795  15.575  1.00 49.84           O  
ATOM     19  N1   rG R   1      95.826 -36.601  14.053  1.00 51.01           N  
ATOM     20  C2   rG R   1      95.261 -37.624  13.321  1.00 52.69           C  
ATOM     21  N2   rG R   1      94.284 -37.269  12.470  1.00 51.75           N  
ATOM     22  N3   rG R   1      95.605 -38.902  13.435  1.00 50.38           N  
01234567890123456789012345678901234567890123456789012345678901234567890123456789	
0         1         2         3         4         5         6         7  

#Convert 1 to 2:
#					ATOM      2  C5*  rG A   1      40.311  50.542  48.455  1.00 33.79           C  (Rosetta)
#					ATOM      2  C5' G   B   1      40.311  50.542  48.455  1.00 33.79           C  (ModeRNA)

#					ATOM     22  O1P  rG A   1      41.587  49.843  50.847  1.00 43.84           O  (Rosetta)
#					ATOM     22  OP1 G   B   1      41.587  49.843  50.847  1.00 43.84           O  (ModeRNA)

#		   		01234567890123456789012345678901234567890123456789012345678901234567890123456789	
#					0         1         2         3         4         5         6         7     
'''
