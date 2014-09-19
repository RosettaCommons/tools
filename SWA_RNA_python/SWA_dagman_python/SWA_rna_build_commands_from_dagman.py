#!/usr/bin/python

import subprocess
from os import listdir
from os.path import exists, isfile, dirname, basename, expanduser, expandvars
from sys import exit, argv
import string

from SWA_rna_build_dag_parse import swa_rna_build_dag_parse

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

#CONDOR_SAMPLER_DIR = 'CONDOR/SAMPLER/'
#CONDOR_CLUSTERER_DIR = 'CONDOR/CLUSTERER/'

#CONDOR_SAMPLER_FILES=[
#    item for item in listdir(CONDOR_SAMPLER_DIR) if isfile(CONDOR_SAMPLER_DIR+item)]
#CONDOR_CLUSTERER_FILES=[
#    item for item in listdir(CONDOR_CLUSTERER_DIR) if isfile(CONDOR_CLUSTERER_DIR+item)]

#CONDOR_FILE_DICT={
#    CONDOR_SAMPLER_DIR : CONDOR_SAMPLER_FILES,
#    CONDOR_CLUSTERER_DIR : CONDOR_CLUSTERER_FILES
#}


#for CONDOR_DIR, CONDOR_FILES in CONDOR_FILE_DICT.iteritems():
#    for condor_filename in CONDOR_FILES:
#        print 'Writing command to ', JOBS_DIR + condor_filename.split('.')[0]
#        condor_file = open( CONDOR_DIR + condor_filename, 'r' )
#        lines = condor_file.readlines()
#        executable = lines[0].split(' = ')[1].strip('\n')
#        arguments = lines[1].split(' = ')[1].strip('\n')
#        condor_file.close()
#
#        command = executable + ' ' + arguments
#        command += ' ' + ' '.join(extra_flags)

#        job_file = open( JOBS_DIR + condor_filename.split('.')[0], 'w' )
#        job_file.write( command )
#        job_file.close()


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


