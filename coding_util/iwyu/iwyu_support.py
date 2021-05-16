#!/usr/bin/env python

'''iwyu_support.py - A library file containing support functions shared between apply_iwyu.py and run_iwyu.py'''

from __future__ import print_function

import sys, os

#These are the clang commandline flags for the cat=test debug mode (You can update them by just copy-pasting from a regular Clang compile).
commandline_flags_debug_linux = '''-c -std=c++11 -isystem external/boost_submod/ -isystem external/ -isystem external/include/ -isystem external/dbio/ -isystem external/libxml2/include -isystem external/rdkit -isystem external/cxxtest/ -pipe -Qunused-arguments -DUNUSUAL_ALLOCATOR_DECLARATION -ftemplate-depth-256 -stdlib=libstdc++ -W -Wall -Wextra -pedantic -Werror -Wno-long-long -Wno-strict-aliasing -O0 -g -Wno-unused-function -DBOOST_ERROR_CODE_HEADER_ONLY -DBOOST_SYSTEM_NO_DEPRECATED -DBOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS -DBOOST_DISABLE_THREADS -DPTR_STD -Iexternal/cxxtest -I. -Isrc -Iexternal/include -Itest -Isrc/platform/linux/64/clang -Isrc/platform/linux/64 -Isrc/platform/linux'''.split()

def get_commandline_flags():
    '''Return the command line flags for the Clang++ run on this system.'''
    # Currently, only support Linux Debug mode.
    # We could add some fancy platform parsing, but that's not needed at the moment.
    if True:
        flags = commandline_flags_debug_linux
    return flags + [
        "-x", "c++",           # I'm not sure why, but this causes clang to (correctly) error out when it otherwise wouldn't
        "-ferror-limit=1",     # Stop on first error
        "-DIWYU_SCAN",         # Special flag to say we're in the IWYU_SCAN environment
        "-D_CXXTEST_HAVE_STD", # A needed define added by the testings system scripts
        "-Wno-unused-variable", # Some headers have unused variables
        ]


# May need to be updated for additional include directories in command line
def check_include_file_exists(filename):
    '''We assume we're running in the main/source directory'''
    return ( os.path.exists( 'src/' + filename ) or
        os.path.exists( 'test/' + filename ) or
        os.path.exists( 'external/include/' + filename ) or
        os.path.exists( 'external/' + filename ) or
        os.path.exists( 'external/boost_submod/' + filename ) or
        os.path.exists( 'external/dbio/' + filename ) or
        os.path.exists( 'external/rdkit/' + filename ) or
        os.path.exists( 'external/libxml2/include/' + filename ) or
        os.path.exists( 'external/cxxtest/' + filename )
    )
