#!/usr/bin/env python

################################################################################

import subprocess
from os import listdir
from os.path import exists, isfile, dirname, basename, expanduser, expandvars
from sys import exit, argv
import string

from SWA_rna_build_dag_parse import swa_rna_build_dag_parse

def replace_args( command ):
    replace_args_dict = {}
    replace_args_dict[ '-sample_virtual_ribose_list' ] = '-sample_virtual_sugar_list'
    replace_args_dict[ 'rna_sample_virtual_ribose' ] = 'rna_sample_virtual_sugar'
    for old_arg, new_arg in replace_args_dict.iteritems():
        if old_arg in command.split():
            command = command.replace( old_arg, new_arg )
    return command



extra_flags = [
    '-sampler_perform_phosphate_pack false',
    '-allow_virtual_side_chains false'
]

JOBS_DIR = 'JOBS/'
SCRIPTS_DIR = 'SCRIPTS/'
SCRIPTS_PRE_DIR = 'SCRIPTS/PRE/'
SCRIPTS_POST_DIR = 'SCRIPTS/POST/'

if not exists( JOBS_DIR ):
    subprocess.call( 'mkdir ' + JOBS_DIR, shell=True )
if not exists( SCRIPTS_DIR ):
    subprocess.call( 'mkdir ' + SCRIPTS_DIR, shell=True )
if not exists( SCRIPTS_PRE_DIR ):
    subprocess.call( 'mkdir ' + SCRIPTS_PRE_DIR, shell=True )
if not exists( SCRIPTS_POST_DIR ):
    subprocess.call( 'mkdir ' + SCRIPTS_POST_DIR, shell=True )

###############################################################################

JOBS, SCRIPTS_PRE, SCRIPTS_POST = swa_rna_build_dag_parse()

for job_name, job_ifname in JOBS.iteritems():
    job_ofname = job_ifname.split('/')[-2] + '_' + job_ifname.split('/')[-1].split('.')[0]
    print 'Writing command to ', JOBS_DIR + job_ofname
    job_src = open( job_ifname, 'r' )
    lines = job_src.readlines()
    executable = lines[0].split(' = ')[1].strip('\n')
    arguments = lines[1].split(' = ')[1].strip('\n')
    job_src.close()

    command = executable + ' ' + arguments
    command += ' ' + ' '.join(extra_flags)

    command = replace_args( command )

    job_fout = open( JOBS_DIR + job_ofname, 'w' )
    job_fout.write( command )
    job_fout.close()

for script_name, script in SCRIPTS_PRE.iteritems():
    print 'Writing command to ', SCRIPTS_PRE_DIR + script_name
    script_fout = open( SCRIPTS_PRE_DIR + script_name, 'w' )
    script_fout.write( script )
    script_fout.close()

for script_name, script in SCRIPTS_POST.iteritems():
    print 'Writing command to ', SCRIPTS_POST_DIR + script_name
    script_fout = open( SCRIPTS_POST_DIR + script_name, 'w' )
    script_fout.write( script )
    script_fout.close()


