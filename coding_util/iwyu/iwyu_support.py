#!/usr/bin/env python

'''iwyu_support.py - A library file containing support functions shared between apply_iwyu.py and run_iwyu.py'''

from __future__ import print_function

import sys, os

#These are the clang commandline flags for debug mode, stripped of warning issues (You can update them by just copy-pasting from a regular Clang compile).
commandline_flags_debug_linux = '''-c -std=c++11 -isystem external/boost_submod/ -isystem external/ -isystem external/include/ -isystem external/dbio/ -isystem external/rdkit -isystem external/libxml2/include -isystem external/cxxtest/ -pipe -Qunused-arguments -DUNUSUAL_ALLOCATOR_DECLARATION -ftemplate-depth-256 -stdlib=libstdc++ -Wno-long-long -Wno-strict-aliasing -O0 -g -DBOOST_ERROR_CODE_HEADER_ONLY -DBOOST_SYSTEM_NO_DEPRECATED -DBOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS -DPTR_STD -Isrc -I./ -Itest/ -Isrc/platform/linux'''.split()


def get_commandline_flags():
    '''Return the command line flags for the Clang++ run on this system.'''
    # Currently, only support Linux Debug mode.
    # We could add some fancy platform parsing, but that's not needed at the moment.
    if True:
        flags = commandline_flags_debug_linux
    # Stop at first error, add the IWYU specific define,
    # Add a define that's injected by the CXXTEST suite generating script
    return flags + ["-ferror-limit=1","-DIWYU_SCAN","-D_CXXTEST_HAVE_STD"]


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
