#!/usr/bin/python

import subprocess
from os import listdir
from os.path import exists, isfile

JOBS_DIR = 'JOBS/'
SCRIPTS_DIR = 'SCRIPTS/'
SCRIPTS_PRE_DIR = 'SCRIPTS/PRE/'
SCRIPTS_POST_DIR = 'SCRIPTS/POST/'


JOB_FILES=[
    item for item in listdir(JOBS_DIR) if isfile(JOBS_DIR+item)]


# JOB 1
JOB = 'REGION_0_1_START_FROM_REGION_0_0'
print "Running command: source " + JOBS_DIR + JOB
subprocess.call( 'source ' + JOBS_DIR + JOB, shell=True )

# JOB 2
JOB = 'REGION_2_0_START_FROM_REGION_0_0'
print "Running command: source " + JOBS_DIR + JOB
subprocess.call( 'source ' + JOBS_DIR + JOB, shell=True )

# JOB 3
JOB = 'REGION_1_0_START_FROM_REGION_0_0'
print "Running command: source " + JOBS_DIR + JOB
subprocess.call( 'source ' + JOBS_DIR + JOB, shell=True )

# JOB 4
JOB = 'REGION_0_2_START_FROM_REGION_0_0'
print "Running command: source " + JOBS_DIR + JOB
subprocess.call( 'source ' + JOBS_DIR + JOB, shell=True )

