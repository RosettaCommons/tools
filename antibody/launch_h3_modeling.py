#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   launch_h3_modeling.py
## @brief  linker script between antibody.py and H3 modeling, for use with multi-template grafting
## @author Nick Marze (nickmarze@gmail.com)

import os, sys, re, json, commands, shutil

from optparse import OptionParser, IndentedHelpFormatter

_script_path_ = os.path.dirname( os.path.realpath(__file__) )

def main(args):
    ''' Script for sending structures grafted by antibody.py as input to Rosetta H3 modeling,
		with functionality for multiple grafting templates
    '''
    parser = OptionParser(usage="usage: %prog [OPTIONS] [TESTS]")
    parser.set_description(main.__doc__)

    parser.add_option('--nstruct-total',
        default=1000, type="int",
        help="Specify total number of decoys to generate in antibody modeling.",
    )

    parser.add_option('--n-grafted-structures',
        default=1, type="int",
        help="Specify number of grafted models to use as input.",
	)

    parser.add_option('--split-ratio',
        action="store", default='1',
        help="Specify decoy split between input grafted models; default is an even split. Use comma-delimited ratio.",
    )

    parser.add_option('--prefix',
        action="store", default='grafting/',
        help="Prefix for output files (directory name). Default is grafting/.",
    )

    parser.add_option('--rosetta-database',
        action="store", default=None,
        help="Specify path of rosetta database dir. Default is '$ROSETTA/main/database', and if not found corresponding steps will be skipped.",
    )

    parser.add_option('--rosetta-bin',
        action="store", default=None,
        help="Specify path to '/bin' dir with antibody_H3. Default is '$ROSETTA/main/source/bin', and if not found corresponding steps will be skipped.",
    )
	
    parser.add_option('--rosetta-platform',
		action="store", default=None,
        help="Specify full extra+compier+build type for rosetta biniaries found in --rosetta-bin. For example use static.linuxgccrelease for static build on Linux. Default is dynamic release build of current OS",
    )

    parser.add_option('--flags-file',
        action="store", default=None,
        help="Specify path of antibody_H3 flags file.",
    )


    print "starting now..."

    #os.getcwd() #os.chdir()
    global script_dir
    script_dir = os.path.dirname(__file__)

    #parse options
    ( options, args ) = parser.parse_args( args = args[1:] )
	
    if not options.flags_file: options.flags_file = script_dir + '/abH3.flags'
    if not options.rosetta_bin:
        if 'ROSETTA' in os.environ:
            options.rosetta_bin = os.path.abspath(os.environ['ROSETTA']) + '/main/source/bin'
        else:
            print 'rosetta binary not found; please specify with --rosetta-database... exiting...'
            sys.exit(1)
    if not options.rosetta_database:
        if 'ROSETTA' in os.environ:
            options.rosetta_database = os.path.abspath(os.environ['ROSETTA']) + '/main/database'
        else:
            print 'rosetta database not found; please specify with --rosetta-bin... exiting...'
            sys.exit(1)
    options.rosetta_bin = os.path.abspath( options.rosetta_bin )
    options.rosetta_database = os.path.abspath( options.rosetta_database )
    if not options.rosetta_platform:
        if sys.platform.startswith("linux"): options.rosetta_platform = 'linuxgccrelease'
        elif sys.platform == "darwin" : options.rosetta_platform = 'macosgccrelease'
        else: options.rosetta_platform = '_unknown_'
    nstruct_total = options.nstruct_total
    n_grafted_structures = options.n_grafted_structures
    split_ratio = options.split_ratio

    #convert decoy ratio to list of nstructs
    split_ratio_list = split_ratio.split( ',' )
    while len( split_ratio_list ) < n_grafted_structures:
        split_ratio_list.append( split_ratio_list[-1] )

    scaling_factor = 0
    for i in range( 0, n_grafted_structures ):
        scaling_factor = scaling_factor + int( split_ratio_list[i] )
    if scaling_factor == 0: scaling_factor = 1

    nstruct_list = [ ( ( int( j ) * nstruct_total ) / scaling_factor ) for j in split_ratio_list ]

    print nstruct_list

    #generate shell scripts
    for k in range( 0, n_grafted_structures ):
        input_file = options.prefix + '/model.' + str(k) + '.pdb'
        constraint_file = options.prefix + '/cter_constraint.' + str( k )
        flags_file = options.flags_file
        database = options.rosetta_database
        model_h3 = options.rosetta_bin + '/antibody_H3.' + options.rosetta_platform
        script = '#!/bin/sh\n' + model_h3 + ' @flags ' + flags_file + ' -database ' + database + ' -s ' + input_file + ' -constraints:cst_file ' + constraint_file + ' -nstruct ' + str( nstruct_list[k] ) + '\n'
        print script
        f=open( options.prefix + 'run_h3_modeling.' + str(k) + '.sh', 'w')
        f.write( script )


if __name__ == "__main__": main(sys.argv)