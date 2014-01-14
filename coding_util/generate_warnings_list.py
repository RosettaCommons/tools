#!/usr/bin/python
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
from os import rename, remove
from subprocess import call
from argparse import ArgumentParser

# Parse arguments.
parser = ArgumentParser(description=__doc__)
parser.add_argument('output_filename', nargs="?", default='warnings.list',
                    help='the output filename for the list of warnings')
parser.add_argument('-c', '--cxx', default='clang', choices=['clang', 'gcc'],
                    help='the compiler to use')
parser.add_argument('-j', type=int, default=1,
                    help='number of processors to use')
parser.add_argument('--mode', default='debug',
                    choices=['debug', 'release'],
                    help='the build mode')
parser.add_argument('--mute', action='store_true',
                    help='flag to mute Scons output during Rosetta build')
parser.add_argument('--settings_directory',
                    default='../../main/source/tools/build',
                    help='directory where the "user.settings" file is located')
parser.add_argument('--source_directory',
                    default='../../main/source',
                    help='directory from which to compile Rosetta')
parser.add_argument('-f', '--format', default='none',
                    choices=['none', 'html', 'wiki'],
                    help='which format for the output list')
args = parser.parse_args()

# Set constants.
COMMAND = ['python', 'scons.py', 'mode=' + args.mode, 'bin', 'cxx=' + args.cxx,
           '-j' + str(args.j)]
RAW_WARNINGS_FILE = open('raw_warnings.tmp', "w")
if args.mute:
    SCREEN_OUTPUT_FILE = open('junk.tmp', "w")
else:
    SCREEN_OUTPUT_FILE = None

# Subroutines
def restore_settings():
    remove(args.settings_directory + '/user.settings')
    rename(args.settings_directory + '/user.settings.temporary',
           args.settings_directory + '/user.settings')

def clean_up():
    pass


# Move user.settings to temporary file.
print 'saving original "user.settings" data...'
rename(args.settings_directory + '/user.settings',
       args.settings_directory + '/user.settings.temporary')

# Generate new user.settings file.
print 'generating new settings:'

settings = 'settings = {"user" : {"prepends" : {}, "appends" : {}, ' + \
    '"overrides" : {"flags" : {"warn" : ['
if args.cxx == 'clang':
    settings += '"Weverything", "fno-caret-diagnostics", ' + \
                '"fno-color-diagnostics", "fno-diagnostics-fixit-info", '
else:  # cxx=gcc
    settings += '"Wall", "Wextra", "pedantic", '
settings += '], }}, "removes" : {}, }}\n'

print settings
with open(args.settings_directory + '/user.settings', "w") as f:
    f.write(settings)

# Delete build directories
print 'cleaning build directories for fresh build...'
#rmdir

# Move to source directory and compile Rosetta with new settings.
print 'building Rosetta with Scons...'
try:
    return_code = call(COMMAND, stdout=SCREEN_OUTPUT_FILE,
                       stderr=RAW_WARNINGS_FILE, cwd=args.source_directory)
    if return_code < 0:
        print >>sys.stderr, 'Scons terminated with the following signal:',
        print >>sys.stderr, -return_code
    else:
        print 'Scons has terminated.'
except OSError as e:
    print >>sys.stderr, 'Compiling of Rosetta failed:', e
except KeyboardInterrupt:
    RAW_WARNINGS_FILE.close()
    if SCREEN_OUTPUT_FILE is not None:
        SCREEN_OUTPUT_FILE.close()
    restore_settings()
    clean_up()
    exit('\nBuild cancelled; original "user.settings" restored.')

# Restore original user.settings file.
print 'restoring original "user.settings" data...'
restore_settings()

# Parse RAW_WARNINGS_FILE.
raw_data = RAW_WARNINGS_FILE.readlines()
print raw_data[-1]  # last element in list
if raw_data[-1].startswith('scons: '):
    print 'Rosetta was NOT successfully compiled. ',
    print 'Re-run with a non-broken version of Rosetta.'
else:
    print 'Rosetta was successfully compiled.'

# Finish.
clean_up()
print 'A list of warnings can be found',
print 'in the "' + args.output_filename + '" file',
print 'in the current directory.'