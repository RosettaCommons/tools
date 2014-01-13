#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available
# (c) under license.
# (c) The Rosetta software is developed by the contributing members of the
# (c) Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington UW
# (c) TechTransfer, email: license@u.washington.edu.

"""Brief:   This script performs a clean compile in the specified compiler and
         generates a list of warnings along with statistics.

Author:  Jason W. Labonte

"""

# Imports
import sys
from argparse import ArgumentParser


# Parse arguments.
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-c', '--cxx', default='clang',
                    help='the compiler to use')
parser.add_argument('-j', type=int, default=1,
                    help='number of processors to use')
args = parser.parse_args()