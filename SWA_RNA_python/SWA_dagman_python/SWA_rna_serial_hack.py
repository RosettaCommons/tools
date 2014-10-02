#!/usr/bin/python

############################################################################

import subprocess
from os import listdir
from os.path import exists, isfile, expandvars

############################################################################

def fill_queue( queue, job=None, script_pre=None, script_post=None ):
    queue[ 'jobs' ].append( job )
    queue[ 'pre_process' ].append( script_pre )
    queue[ 'post_process' ].append( script_post )
    return queue

def set_process_id( command_src, process_id ):
    command_fin = open( command_src, 'r' )
    command = ' '.join( command_fin.readlines() )
    command_fin.close()
    if '$(Process)' in command:
        command_fout = open( command_src, 'w' )
        command_fout.write( command.replace( '$(Process)', '%i' % process_id ) ) 
        command_fout.close()
    return True

def filter_cluster_out( jobname ):
    if jobname == 'CLUSTERER_REGION_0_1_cluster':
        cluster_out_fname = 'region_0_1_sample.cluster.out'

    elif jobname == 'CLUSTERER_REGION_2_0_cluster':
        cluster_out_fname = 'region_2_0_sample.cluster.out'

    else: return


    if exists( cluster_out_fname ):
        print "cluster_out_fname: "+cluster_out_fname+" exists!!!"
        fin = open( cluster_out_fname, 'r' )
        lines = fin.readlines()
        for line in lines:
            if 'ANNOTATED_SEQUENCE:' in line.split()[0]:
                tag = line.split()[-1]
                break
        fin.close()
        fout = open( cluster_out_fname, 'w')
        for line in lines:
            if 'S_' in line.split()[-1] and tag not in line.split()[-1]: break
            fout.write(line)
        fout.close()

        return



#############################################################################

build_commands = '$ROSETTA/tools/SWA_RNA_python/SWA_dagman_python/SWA_rna_build_commands_from_dagman.py'
subprocess.call( 'python ' + build_commands , shell=True )

JOBS_DIR = 'JOBS/'
SCRIPTS_DIR = 'SCRIPTS/'
SCRIPTS_PRE_DIR = 'SCRIPTS/PRE/'
SCRIPTS_POST_DIR = 'SCRIPTS/POST/'


JOB_FILES=[
    item for item in listdir(JOBS_DIR) if isfile(JOBS_DIR+item)]

queue={
   'jobs'        : [],
   'pre_process' : [],
   'post_process': []
}

#############################################################################
### SAMPLING ROUND 1
#############################################################################

# SAMPLE 0 1
JOB = 'SAMPLER_REGION_0_1_START_FROM_REGION_0_0'
SCRIPT_POST = 'REGION_0_1_START_FROM_REGION_0_0'
queue = fill_queue( queue, job=JOB, script_post=SCRIPT_POST )

# SAMPLE 2 0
JOB = 'SAMPLER_REGION_2_0_START_FROM_REGION_0_0'
SCRIPT_POST = 'REGION_2_0_START_FROM_REGION_0_0'
queue = fill_queue( queue, job=JOB, script_post=SCRIPT_POST )

# SAMPLE 1 0
JOB = 'SAMPLER_REGION_1_0_START_FROM_REGION_0_0'
SCRIPT_POST = 'REGION_1_0_START_FROM_REGION_0_0'
queue = fill_queue( queue, job=JOB, script_post=SCRIPT_POST )

# SAMPLE 0 2
JOB = 'SAMPLER_REGION_0_2_START_FROM_REGION_0_0'
SCRIPT_POST = 'REGION_0_2_START_FROM_REGION_0_0'
queue = fill_queue( queue, job=JOB, script_post=SCRIPT_POST )


############################################################################
### CLUSTERING ROUND 1
############################################################################

# CLUSTER 0 1
JOB ='CLUSTERER_REGION_0_1_cluster'
queue = fill_queue( queue, job=JOB )

# CLUSTER 2 0
JOB ='CLUSTERER_REGION_2_0_cluster'
queue = fill_queue( queue, job=JOB )


############################################################################                                                                 ### VIRT SAMPLING ROUND 1
############################################################################

# VIRT 0 1 FOR 0 2 FROM 0 1
JOB =         'VIRT_RIBOSE_SAMPLER_REGION_0_1_FOR_REGION_0_2_START_FROM_REGION_0_1'
SCRIPT_PRE =  'VIRT_RIBOSE_SAMPLER_REGION_0_1_FOR_REGION_0_2_START_FROM_REGION_0_1'
SCRIPT_POST = 'VIRT_RIBOSE_SAMPLER_REGION_0_1_FOR_REGION_0_2_START_FROM_REGION_0_1'
queue = fill_queue( queue, job=JOB, script_pre=SCRIPT_PRE, script_post=SCRIPT_POST )

# VIRT 2 0 FOR 1 0 FROM 2 0
JOB =         'VIRT_RIBOSE_SAMPLER_REGION_2_0_FOR_REGION_1_0_START_FROM_REGION_2_0'
SCRIPT_PRE =  'VIRT_RIBOSE_SAMPLER_REGION_2_0_FOR_REGION_1_0_START_FROM_REGION_2_0'
SCRIPT_POST = 'VIRT_RIBOSE_SAMPLER_REGION_2_0_FOR_REGION_1_0_START_FROM_REGION_2_0'
queue = fill_queue( queue, job=JOB, script_pre=SCRIPT_PRE, script_post=SCRIPT_POST )

# VIRT 0 1 FOR 2 1 FROM 0 1 AND 0 0
JOB =         'VIRT_RIBOSE_SAMPLER_REGION_0_1_FOR_REGION_2_1_START_FROM_REGION_0_1_AND_0_0'
SCRIPT_PRE =  'VIRT_RIBOSE_SAMPLER_REGION_0_1_FOR_REGION_2_1_START_FROM_REGION_0_1_AND_0_0'
SCRIPT_POST = 'VIRT_RIBOSE_SAMPLER_REGION_0_1_FOR_REGION_2_1_START_FROM_REGION_0_1_AND_0_0'
queue = fill_queue( queue, job=JOB, script_pre=SCRIPT_PRE, script_post=SCRIPT_POST )

# VIRT 2 0 FOR 2 1 FROM 0 0 AND 2 0
JOB =         'VIRT_RIBOSE_SAMPLER_REGION_2_0_FOR_REGION_2_1_START_FROM_REGION_0_0_AND_2_0'
SCRIPT_PRE =  'VIRT_RIBOSE_SAMPLER_REGION_2_0_FOR_REGION_2_1_START_FROM_REGION_0_0_AND_2_0'
SCRIPT_POST = 'VIRT_RIBOSE_SAMPLER_REGION_2_0_FOR_REGION_2_1_START_FROM_REGION_0_0_AND_2_0'
queue = fill_queue( queue, job=JOB, script_pre=SCRIPT_PRE, script_post=SCRIPT_POST )


############################################################################                                                                 ### SAMPLING ROUND 2
############################################################################

# SAMPLE 0 2 FROM 0 1
JOB = 'SAMPLER_REGION_0_2_START_FROM_REGION_0_1'
SCRIPT_PRE =  'REGION_0_2_START_FROM_REGION_0_1'
SCRIPT_POST = 'REGION_0_2_START_FROM_REGION_0_1'
queue = fill_queue( queue, job=JOB, script_pre=SCRIPT_PRE, script_post=SCRIPT_POST )

# SAMPLE 1 0 FROM 2 0
JOB = 'SAMPLER_REGION_1_0_START_FROM_REGION_2_0'
SCRIPT_PRE =  'REGION_1_0_START_FROM_REGION_2_0'
SCRIPT_POST = 'REGION_1_0_START_FROM_REGION_2_0'
queue = fill_queue( queue, job=JOB, script_pre=SCRIPT_PRE, script_post=SCRIPT_POST )

# SAMPLE 2 1 FROM 0 1 AND 0 0
JOB = 'SAMPLER_REGION_2_1_START_FROM_REGION_0_1_AND_0_0'
SCRIPT_PRE =  'REGION_2_1_START_FROM_REGION_0_1_AND_0_0'
SCRIPT_POST = 'REGION_2_1_START_FROM_REGION_0_1_AND_0_0'
queue = fill_queue( queue, job=JOB, script_pre=SCRIPT_PRE, script_post=SCRIPT_POST )

# SAMPLE 2 1 FROM 0 0 AND 2 0
JOB = 'SAMPLER_REGION_2_1_START_FROM_REGION_0_0_AND_2_0'
SCRIPT_PRE =  'REGION_2_1_START_FROM_REGION_0_0_AND_2_0'
SCRIPT_POST = 'REGION_2_1_START_FROM_REGION_0_0_AND_2_0'
queue = fill_queue( queue, job=JOB, script_pre=SCRIPT_PRE, script_post=SCRIPT_POST )


############################################################################                                                                 ### CLUSTERING ROUND 2
############################################################################

# CLUSTER 0 2
JOB ='CLUSTERER_REGION_0_2_cluster'
queue = fill_queue( queue, job=JOB )

# CLUSTER 1 0
JOB ='CLUSTERER_REGION_1_0_cluster'
queue = fill_queue( queue, job=JOB )

# CLUSTER 2 1
JOB ='CLUSTERER_REGION_2_1_cluster'
queue = fill_queue( queue, job=JOB )


############################################################################                                                                 ### CLUSTERING FINAL
############################################################################
# CLUSTER FINAL
JOB ='CLUSTERER_REGION_FINAL_cluster'
queue = fill_queue( queue, job=JOB )


############################################################################
### RUN QUEUED COMMANDS IN SERIAL
############################################################################

for ii in xrange( len( queue[ 'jobs' ] ) ):

    ### PRE PROCESSING SCRIPT
    if queue[ 'pre_process' ][ ii ]:
        command_src = SCRIPTS_PRE_DIR + queue[ 'pre_process' ][ ii ]
        print "Running command: source " + command_src
        subprocess.call( 'source ' + command_src, shell=True )

    ### JOB
    if queue[ 'jobs' ][ ii ]:
        command_src = JOBS_DIR + queue[ 'jobs' ][ ii ]
        dummy = 0
        set_process_id( command_src, dummy )
        print "Running command: " + command_src
        subprocess.call( 'source ' + command_src, shell=True )
        if( queue[ 'jobs' ][ ii ] == 'CLUSTERER_REGION_0_1_cluster' or
            queue[ 'jobs' ][ ii ] == 'CLUSTERER_REGION_2_0_cluster' ):
            filter_cluster_out( queue[ 'jobs' ][ ii ])


    ### POST PROCESSING SCRIPT
    if queue[ 'post_process' ][ ii ]:
        command_src = SCRIPTS_POST_DIR + queue[ 'post_process' ][ ii ]
        print "Running command: source " + command_src
        subprocess.call( 'source ' + command_src, shell=True )










