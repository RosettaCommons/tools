#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


"""Brief:   This script performs a clean compile in the specified compiler and
         generates a list of warnings along with statistics.

Details: This script operates by creating a user.settings file intended to
         have Scons generate all warnings.  It performs a clean build of
         Rosetta and then parses the generated list of warnings on the
         assumption that all warnings output lines end in a square-bracketed
         designation of the warning type.  Statistics are generated and results
         are output in the requested format.  The script can only compile
         Rosetta using clang or gcc.

Note:    Sometimes Scons will return corrupted/nonsensical warning output, so
         it is impossible to get a perfect count.

Params:  ./generate_warnings_list.py -h displays all parameter options.

Example: ./generate_warnings_list.py

Author:  Jason W. Labonte

"""

# Imports
import sys
from os import rename, remove
from subprocess import call
from shutil import rmtree
from argparse import ArgumentParser


# Parse arguments.
parser = ArgumentParser(description=__doc__)
parser.add_argument('output_filename', nargs="?", default='warnings',
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
parser.add_argument('--keep_raw', action='store_true',
                    help='do not delete the raw data file of warnings output')
parser.add_argument('--parse_only', action='store_true',
                    help='generate list and counts from already-existing raw' +
                    'data file')
args = parser.parse_args()
if args.parse_only:
    args.keep_raw = True


# Subroutines
def restore_settings():
    remove(args.settings_directory + '/user.settings')
    rename(args.settings_directory + '/user.settings.temporary',
           args.settings_directory + '/user.settings')


def clean_up():
    try:
        remove('junk.tmp')
    except OSError:
        pass
    if not args.keep_raw:
        try:
            remove('raw_warnings.tmp')
        except OSError:
            pass


def build_rosetta(out, err):
    print 'building Rosetta with Scons...'
    command = ['python', 'scons.py', 'mode=' + args.mode, 'bin',
               'cxx=' + args.cxx, '-j' + str(args.j)]
    try:
        return_code = call(command, stdout=out, stderr=err,
                           cwd=args.source_directory)
        if return_code < 0:
            print >>sys.stderr, 'Scons terminated with the following signal:',
            print >>sys.stderr, -return_code
        else:
            print 'Scons has terminated.'
    except OSError as e:
        print >>sys.stderr, 'Compiling of Rosetta failed:', e
    except KeyboardInterrupt:
        restore_settings()
        clean_up()
        exit('\nBuild canceled; original "user.settings" restored.')


def parse_gcc_warnings(raw_data):
    pass


def parse_clang_warnings(raw_data):
    pass


# Main
if not args.parse_only:
    # Move user.settings to temporary file.
    print 'saving original "user.settings" data...'
    rename(args.settings_directory + '/user.settings',
           args.settings_directory + '/user.settings.temporary')

    # Generate new user.settings file.
    print 'generating new settings:'

    settings = 'settings = {"user" : {"prepends" : {}, "appends" : ' +\
               '{"flags" : {"warn" : ['
    if args.cxx == 'clang':
        settings += '"Wall", "fno-caret-diagnostics", ' + \
                    '"fno-color-diagnostics", "fno-diagnostics-fixit-info",'
    else:  # cxx=gcc
        settings += '"Wall", "Wextra", "pedantic", ' + \
                    '"fdiagnostics-show-option", "fmessage-length=0"'
    settings += '], }}, "overrides" : {}, "removes" : {}, }}\n'

    print settings
    with open(args.settings_directory + '/user.settings', "w") as f:
        f.write(settings)

    # Delete build directories
    print 'cleaning build directories for fresh build...'
    try:
        rmtree(args.source_directory + '/build/external')
    except OSError:
        pass
    try:
        rmtree(args.source_directory + '/build/src')
    except OSError:
        pass

    # Move to source directory and compile Rosetta with new settings.
    with open('raw_warnings.tmp', "w") as raw_warnings_file:
        if args.mute:
            with open('junk.tmp', "w") as screen_output_file:
                build_rosetta(screen_output_file, raw_warnings_file)
        else:
            build_rosetta(None, raw_warnings_file)

    # Restore original user.settings file.
    print 'restoring original "user.settings" data...'
    restore_settings()


# Parse raw_warnings file.
with open('raw_warnings.tmp', "r") as f:
    raw_data = f.readlines()
if raw_data[-1].startswith('scons: '):  # last element in list
    print 'Rosetta was NOT successfully compiled. ',
    print 'Re-run with a non-broken version of Rosetta.'
else:
    print 'Rosetta was successfully compiled. ',
    print 'parsing raw data...'

data = [line for line in raw_data if line.endswith(']\n')]
#        and line.startswith('src/') and not line.startswith('src/ObjexxFCL/')]
sorted_data = {}
bad_lines = 0
for line in data:
    try:
        # Get string between "[" and "]".
        warning_type = line[line.rindex('[-') + 1 : -2]
    except ValueError:
        bad_lines += 1
        continue
    if warning_type in sorted_data:
        sorted_data[warning_type].append(line)
    else:
        sorted_data[warning_type] = [line]
print "Note: Scons produced", bad_lines,
print "bad lines that may also have been warnings;",
print "these will not be counted."


# Generate output files.
total_warnings = str(len(data))
with open(args.output_filename + '.list', "w") as list_file:
    with open(args.output_filename + '.counts', "w") as counts_file:
        list_file.write('TOTAL WARNINGS: ' + total_warnings + "\n")
        counts_file.write('TOTAL\t' + total_warnings + "\n")
        for key, warnings in sorted_data.iteritems():
            # List file
            n_warnings = str(len(warnings))
            list_file.write("\n")
            if args.format == 'html':
                list_file.write('<h3>')
            elif args.format == 'wiki':
                list_file.write('===')
            list_file.write(key)
            list_file.write(' (' + n_warnings + "/" + total_warnings + ")")
            if args.format == 'html':
                list_file.write('</h3>\n')
            elif args.format == 'wiki':
                list_file.write('===\n')
            else:
                list_file.write("\n")
            for warning in warnings:
                if args.format == 'html':
                    warning = '<li>' + warning
                elif args.format == 'wiki':
                    warning = '* ' + warning
                list_file.write(warning)

            # Counts file
            counts_file.write(key)
            counts_file.write("\t")
            counts_file.write(n_warnings)
            counts_file.write("\n")

# Finish.
clean_up()
print 'A list of warnings can be found',
print 'in the "' + args.output_filename + '.list" file',
print 'in the current directory.'
print 'A tab-delimited summary of warning counts can be found',
print 'in the "' + args.output_filename + '.summary" file',
print 'in the current directory.'
