#!/usr/bin/env python
from __future__ import print_function
from sys import argv,exit
from os import system,getcwd,popen,devnull
from os.path import basename,dirname,expanduser,exists,expandvars
import string

def Help():
    print argv[0]+' <text file with rosetta command> <outdir> <# jobs>  [# hours]'
    print '  give outdir as 0 and # jobs as 1 to not create separate outdirs.'
    exit()

if len( argv ) < 2:
    Help()

infile = argv[1]

try:
    outdir = argv[2]
except:
    outdir = '0'

try:
    n_jobs = int( argv[3] )
except:
    n_jobs = 1
    if outdir != '0': print 'NEED TO SUPPLY NUMBER OF JOBS'

DO_MPI = False # May need to reactivate for XSEDE.

hostname = ''
hostname_tag = popen( 'hostname' ).readlines()[0]
if hostname_tag.find( 'ls4' ) > -1: hostname = 'lonestar'
if hostname_tag.find( 'DasLab' ) > -1 or hostname_tag.find( 'das' ) > -1: hostname = 'ade'
if hostname_tag.find( 'stampede' ) > -1: hostname = 'stampede'
if hostname_tag.find( 'comet' ) > -1: hostname = 'comet'
if hostname_tag.find( 'sherlock' ) > -1 or hostname_tag.find( 'sh-' ) > -1: hostname = 'sherlock'
if hostname_tag.find( 'biox3' ) > -1: hostname = 'biox3'

queue_cmd = 'qsub'
if hostname in ["stampede", "sherlock", "comet"]:
    queue_cmd = 'sbatch'

if hostname == 'lonestar':
    DO_MPI = True
    tasks_per_node_MPI = 12 # lonestar
if hostname == 'stampede':
    DO_MPI = True
    tasks_per_node_MPI = 16
    account = 'TG-MCB120152'
if hostname == 'comet':
    DO_MPI = True
    tasks_per_node_MPI = 24
    account = expandvars("$COMET_ACCOUNT")
    if not len(account):
        account = 'TG-MCB120152'
if hostname == 'sherlock':
    DO_MPI = True
    tasks_per_node_MPI = n_jobs if n_jobs < 16 else 16
    account = None

save_logs = False
if argv.count( '-save_logs' )>0:
    save_logs = True
    pos = argv.index( '-save_logs' )
    del( argv[pos] )
development = False
if argv.count( '-development' )>0:
    development = True
    pos = argv.index( '-development' )
    del( argv[pos] )

# set name of queue to submit jobs to
queue = 'normal'
user_input_queue = False
if '-queue' in argv:
    idx = argv.index('-queue')
    argv.pop(idx)
    queue = argv.pop(idx)
    user_input_queue = True
elif development is True:
    queue = 'development'
elif hostname in ['comet']:
    queue = 'compute'


nhours = 16
if len( argv ) > 4:
    nhours = int( argv[4] )
    if ( nhours > 168 ):  Help()
    if hostname in ['sherlock', 'comet']:
        nhours = min(nhours, 48)

if not exists( infile ):
    print('Could not find: ', infile)
    exit( 0 )

lines = open(infile).readlines()

# legacy queue systems -- single file handles all submissions.
bsub_file = 'bsubMINI'
condor_file = 'condorMINI'
qsub_file = 'qsubMINI'
sbatch_file = 'sbatchMINI'

if queue_cmd == 'qsub' or queue_cmd == 'sbatch':
    bsub_file = devnull
    condor_file = devnull
if queue_cmd == 'qsub': sbatch_file = devnull
if queue_cmd == 'sbatch': qsub_file = devnull

# MPI based systems require separate submission files for each command.
queue_file_MPI = queue_cmd + 'MPI'
queue_file_MPI_ONEBATCH = queue_cmd+'MPI_ONEBATCH' # testing on stampede?
job_file_MPI_ONEBATCH = 'MPI_ONEBATCH.job'
all_commands_file = 'all_commands.sh'
if hostname in ["stampede", "sherlock", "comet"]:
    queue_file_MPI_ONEBATCH = "job.batch"
if hostname == 'sherlock':
    queue_file_MPI_ONEBATCH = devnull
    job_file_MPI_ONEBATCH = devnull

fid = open( bsub_file,'w')
fid_condor = open( condor_file,'w')
fid_qsub = open( qsub_file,'w')
fid_sbatch = open( sbatch_file, 'w')
fid_all_commands = open( all_commands_file, 'w' )

if DO_MPI:
    fid_queue_MPI = open( queue_file_MPI,'w')
    fid_job_MPI_ONEBATCH  = open( job_file_MPI_ONEBATCH,'w')
    fid_queue_MPI_ONEBATCH = open( queue_file_MPI_ONEBATCH,'w')

tot_jobs = 0

universe = 'vanilla';
fid_condor.write('+TGProject = "TG-MCB090153"\n')
fid_condor.write('universe = %s\n' % universe)
fid_condor.write('notification = never\n')

HOMEDIR = expanduser('~')
CWD = getcwd()

PATH_LIST = expandvars('$PATH').split(':')
QSUB_SUBMIT_CMD = ('qsub.sh' if sum([exists('%s/qsub.sh'%p) for p in PATH_LIST]) else 'qsub')
SBATCH_SUBMIT_CMD = 'sbatch'

qsub_file_dir = 'qsub_files/'
if queue_cmd == 'qsub' and not exists( qsub_file_dir ): system( 'mkdir '+qsub_file_dir )

if queue_cmd == 'sbatch':
    sbatch_file_dir = 'sbatch_files/'
    if not exists( sbatch_file_dir ): system( 'mkdir '+sbatch_file_dir )

if DO_MPI:
    queue_file_dir_MPI = queue_cmd + '_files_MPI/'
    if not exists( queue_file_dir_MPI ): system( 'mkdir '+queue_file_dir_MPI )

command_lines_explicit = []

if save_logs:
    outfile_general = '$(Process).log'
    errfile_general = '$(Process).err'
else:
    outfile_general = '/dev/null'
    errfile_general = '/dev/null'

for line in lines:

    if len(line) == 0: continue
    if line[0] == '#': continue
    #if line[0].split() == []: continue
    command_line = line[:-1]

    cols = command_line.split()
    for i in range( len( cols ) ):
        if cols[i][0] == '@':
            flag_file = cols[i][1:]
            flags = open( flag_file ).readlines()
            new_flags = ''
            for flag in flags:
                if len(flag)>0 and flag[0] != '#':
                    new_flags += ' ' + flag.replace( '\n', '')
            cols[i] = new_flags
            command_line = ' '.join(cols)

    dir = outdir + '/$(Process)/'
<<<<<<< HEAD
    make_outdirs = False
    if outdir == '0':
        assert( n_jobs == 1 )
    elif command_line.find( '-csa_bank_size' ) > -1:
        print "Detected CSA mode"
=======
    if command_line.find( '-csa_bank_size' ) > -1:
        print("Detected CSA mode")
>>>>>>> master
    else:
        command_line = command_line.replace( 'out:file:silent  ','out:file:silent ').replace( '-out:file:silent ', '-out:file:silent '+dir)
        command_line = command_line.replace( '-out::file::silent ', '-out::file::silent '+dir)
        command_line = command_line.replace( '-silent ', '-out:file:silent '+dir)
        command_line = command_line.replace( '-out:file:o ', '-out:file:o '+dir)
        command_line = command_line.replace( '-output_histogram_file ', '-output_histogram_file '+dir)
        command_line = command_line.replace( '-o ', '-o '+dir)
        make_outdirs = True
    #command_line = command_line.replace( '-seed_offset 0', '-seed_offset $(Process)')
    command_line = command_line.replace( '-constant_seed', '-constant_seed -jran $(Process)')
    command_line = command_line.replace( 'macosgcc', 'linuxgcc')
    command_line = command_line.replace( 'Users', 'home')
    command_line = command_line.replace( '~/', HOMEDIR+'/')
    command_line = command_line.replace( '/home/rhiju',HOMEDIR)

    cols = command_line.split()
    if len( cols ) == 0: continue

    EXE = cols[ 0 ]
    if not exists( EXE ):
        rosetta_folder = expandvars("$ROSETTA")
        EXE = rosetta_folder + '/main/source/bin/'+basename(EXE)
    if not exists( EXE ):
        EXE += ".linuxiccrelease"
    if not exists( EXE ):
        EXE += ".linuxclangrelease"
    if not exists( EXE ):
        EXE += ".macosclangrelease"
    if not exists( EXE ):
        EXE = EXE.replace( '/home1/','/work1/')
    assert( exists( EXE ) )

    if not exists( cols[0] ):
        cols[0] = EXE
        command_line = ' '.join(cols)


    if '-total_jobs' in cols:
        pos = cols.index( '-total_jobs' )
        cols[ pos+1 ] = '%d' % n_jobs
        command_line = ' '.join(cols)
    if '-job_number' in cols:
        pos = cols.index( '-job_number' )
        cols[ pos+1 ] = '$(Process)'
        command_line = ' '.join(cols)


    for i in range( n_jobs ):
        if make_outdirs:
            dir_actual = dir.replace( '$(Process)', '%d' % i)
            system( 'mkdir -p '+ dirname(dir_actual) )
        dir_actual = './'

        outfile = outfile_general.replace( '$(Process)', '%d' % i )
        errfile = errfile_general.replace( '$(Process)', '%d' % i )

        command =  'bsub -W %d:0 -o %s -e %s ' % (nhours, outfile, errfile )
        command_line_explicit = command_line.replace( '$(Process)', '%d' % i )
        command += command_line_explicit
        fid.write( command + '\n')

        if queue_cmd == 'qsub':
            # qsub
            pbs_outfile = '/dev/null'
            pbs_errfile = '/dev/null'
            qsub_submit_file = '%s/qsub%d.sh' % (qsub_file_dir, tot_jobs )
            fid_qsub_submit_file = open( qsub_submit_file, 'w' )
            fid_qsub_submit_file.write( '#!/bin/bash\n'  )
            fid_qsub_submit_file.write('#PBS -N %s\n' %  (CWD+'/'+dir_actual[:-1]).replace( '/', '_' ) )
            fid_qsub_submit_file.write('#PBS -o %s\n' % pbs_outfile)
            fid_qsub_submit_file.write('#PBS -e %s\n' % pbs_errfile)
            fid_qsub_submit_file.write('#PBS -m n\n') # no mail
            fid_qsub_submit_file.write('#PBS -M nobody@stanford.edu\n') # no mail
            #fid_qsub_submit_file.write('#PBS -l mem=500Mb\n'  )
            fid_qsub_submit_file.write('#PBS -l walltime=%d:00:00\n\n' % nhours )
            fid_qsub_submit_file.write( 'cd %s\n\n' % CWD )
            fid_qsub_submit_file.write( '%s > %s 2> %s \n' % (command_line_explicit,outfile,errfile) )
            fid_qsub_submit_file.close()

            fid_qsub.write( '%s %s\n' % ( QSUB_SUBMIT_CMD, qsub_submit_file ) )

        if queue_cmd == 'sbatch':

            # sbatch (no mpi)
            queue2 = queue
            if not user_input_queue:
                if hostname in ['comet']:    queue2 = 'shared'
                if hostname in ['sherlock']: queue2 = 'owners'

            job_name = (basename(CWD)+'/'+dir_actual[:-1]).replace( '/', '_' )

            sbatch_submit_file = '%s/job%d.sbatch' % (sbatch_file_dir, tot_jobs )
            fid_sbatch_submit_file = open( sbatch_submit_file, 'w' )
            fid_sbatch_submit_file.write( '#!/bin/bash\n'  )
            fid_sbatch_submit_file.write( '#SBATCH -J %s\n' % job_name )
            fid_sbatch_submit_file.write( '#SBATCH -o %s\n' % outfile )
            fid_sbatch_submit_file.write( '#SBATCH -e %s\n' % errfile )
            fid_sbatch_submit_file.write( '#SBATCH -p %s\n' % queue2 )
            fid_sbatch_submit_file.write( '#SBATCH -t %d:00:00\n' % nhours )
            fid_sbatch_submit_file.write( '#SBATCH -n %d\n' % 1 )
            fid_sbatch_submit_file.write( '#SBATCH -N %d\n' % 1 )
            if account: fid_sbatch_submit_file.write( '#SBATCH -A %s\n' % account )
            fid_sbatch_submit_file.write( 'cd %s\n\n' % CWD )
            fid_sbatch_submit_file.write( '%s\n' % (command_line_explicit) )
            fid_sbatch_submit_file.close()
            fid_sbatch.write( '%s %s\n' % ( SBATCH_SUBMIT_CMD, sbatch_submit_file ) )


        # MPI job file
        if DO_MPI:
            if queue_cmd == 'sbatch':
                fid_job_MPI_ONEBATCH.write( '%s\t%s \n' % (CWD, command_line_explicit ) )
            else:
                fid_job_MPI_ONEBATCH.write( '%s ;;; %s\n' % (CWD, command_line_explicit) )

        # all command lines in one .sh file for bash.
        fid_all_commands.write( 'nohup %s > %s 2> %s &\n' % ( command_line_explicit, outfile, errfile ) )

        command_lines_explicit.append( command_line_explicit )
        tot_jobs += 1

    arguments = ' '.join(cols[1:])

    fid_condor.write('\nexecutable = %s\n' % EXE )
    fid_condor.write('arguments = %s\n' % arguments)
    if save_logs:
        fid_condor.write( 'output = %s\n' % outfile_general )
        fid_condor.write( 'error  = %s\n' % errfile_general )
    fid_condor.write('Queue %d\n' % n_jobs )

if DO_MPI:
    N_MPIJOBS_ONEBATCH =  tot_jobs
    tot_nodes = ( tot_jobs / tasks_per_node_MPI )

    if ( N_MPIJOBS_ONEBATCH % tasks_per_node_MPI != 0 ):  # checking modulo 12 (or whatever the number of cores/node)
        N_MPIJOBS_ONEBATCH += ( tot_jobs/ tasks_per_node_MPI + 1) * tasks_per_node_MPI
        tot_nodes += 1

    count = 0
    for n in range( tot_nodes ):
        # big pain in the ass -- the way we submit jobs via ppserver.py  is
        # failing unless we split the job processor by processor!
        job_submit_file_MPI = '%s/%sMPI%d.job' % (queue_file_dir_MPI, queue_cmd, n )

        fid_job_submit_file_MPI  = open( job_submit_file_MPI, 'w' )
        for m in range( tasks_per_node_MPI ):
            count = count + 1
            if ( count <= tot_jobs ):
                outfile = outfile_general.replace( '$(Process)', '%d' % (count-1) )
                errfile = errfile_general.replace( '$(Process)', '%d' % (count-1) )
                command_line_explicit = command_lines_explicit[ count-1 ]
                if hostname in ["stampede", "sherlock", "comet"]:
                    fid_job_submit_file_MPI.write( '%s\t%s \n' % (CWD, command_line_explicit ) )
                else:
                    command_line_explit +=  ' > %s 2> %s' % (outfile, errfile)
                    fid_job_submit_file_MPI.write( '%s ;;; %s\n' % (CWD, command_line_explicit) )
        fid_job_submit_file_MPI.close()

        if queue_cmd == 'sbatch':
            # sbatch MPI
            jobname= (CWD + '/' + outdir).replace( '/', '_' )
            jobname = jobname[-30:]
            queue_submit_file_MPI = '%s/job%d.batch' % (queue_file_dir_MPI, n )
            fid_queue_submit_file_MPI = open( queue_submit_file_MPI, 'w' )
            fid_queue_submit_file_MPI.write( '#!/bin/bash\n' )
            job_name = (basename(CWD)).replace( '/', '_' )

            # logic from calebgeniess -- may be deprecated
            #queue = 'normal'
            #if hostname in ['comet']:
            #    queue = 'compute'
            #if development:
            #    queue = 'development'
            if hostname in ['sherlock']:
                queue='owners'
            fid_queue_submit_file_MPI.write( '#SBATCH -J %s\n' % job_name )
            fid_queue_submit_file_MPI.write( '#SBATCH -o %s.o%%j\n' % job_name )
            fid_queue_submit_file_MPI.write( '#SBATCH -p %s\n' % queue)
            if development:
                fid_queue_submit_file_MPI.write( '#SBATCH -t 00:10:00\n' )
            else:
                fid_queue_submit_file_MPI.write( '#SBATCH -t %d:00:00\n' % nhours )
            #fid_queue_submit_file_MPI.write( '#SBATCH --mail-user=rhiju@stanford.edu\n' )
            #fid_queue_submit_file_MPI.write( '#SBATCH --mail-type=ALL\n' )
            fid_queue_submit_file_MPI.write( '#SBATCH -n %d\n' % tasks_per_node_MPI )
            fid_queue_submit_file_MPI.write( '#SBATCH -N %d\n' % 1 )
            if account: fid_queue_submit_file_MPI.write( '#SBATCH -A %s\n' % account )
            #fid_queue_submit_file_MPI.write( 'echo $SLURM_NODELIST > nodefile.txt\n' )
            fid_queue_submit_file_MPI.write( 'pp_jobsub.py %s -cluster_name %s -nodelist $SLURM_NODELIST -job_cpus_per_node $SLURM_JOB_CPUS_PER_NODE\n'
                                            % (job_submit_file_MPI, hostname) )
            fid_queue_submit_file_MPI.close()

            fid_queue_MPI.write( 'sbatch %s\n' % queue_submit_file_MPI )
        else:
            # qsub MPI
            jobname= (CWD + '/' + outdir).replace( '/', '_' )
            jobname = jobname[-30:]
            queue_submit_file_MPI = '%s/qsubMPI%d.sh' % (queue_file_dir_MPI, n )
            fid_queue_submit_file_MPI = open( queue_submit_file_MPI, 'w' )
            fid_queue_submit_file_MPI.write( '#!/bin/bash     \n')
            fid_queue_submit_file_MPI.write( '#$ -V  #Inherit the submission environment\n')
            fid_queue_submit_file_MPI.write( '#$ -cwd    # Start job in submission directory\n')
            fid_queue_submit_file_MPI.write( '#$ -N %s   # Job Name\n' % jobname )
            fid_queue_submit_file_MPI.write( '#$ -j y    # Combine stderr and stdout\n')
            fid_queue_submit_file_MPI.write( '#$ -o $JOB_NAME.o$JOB_ID   # Name of the output file\n')
            fid_queue_submit_file_MPI.write( '#$ -pe %dway %d    # Requests X (=12) tasks/node, Y (=12) cores total (Y must be multiples of 12, set X to 12 for lonestar)\n' % (tasks_per_node_MPI, tasks_per_node_MPI) )
            fid_queue_submit_file_MPI.write( '#$ -q %s   # Queue name \n' % queue )
            if nhours == 0: # for testing
                fid_queue_submit_file_MPI.write( '#$ -l h_rt=00:01:00    # Run time (hh:mm:ss)\n' )
            else:
                fid_queue_submit_file_MPI.write( '#$ -l h_rt=%02d:00:00  # Run time (hh:mm:ss)\n' % nhours)
            #fid_queue_submit_file_MPI.write( '#$ -M rhiju@stanford.edu  # Address for email notification\n')
            fid_queue_submit_file_MPI.write( '#$ -m be   # Email at Begin and End of job\n')
            fid_queue_submit_file_MPI.write( 'set -x     # Echo commands, use set echo with csh\n')
            fid_queue_submit_file_MPI.write( 'ibrun mpi_simple_job_submit.py %s  # Run the MPI python\n' % job_submit_file_MPI)
            fid_queue_submit_file_MPI.close()

            fid_queue_MPI.write( 'qsub %s\n' % queue_submit_file_MPI )

    # single batch. does not appear to work...
    # from calebgeniesse -- may need to revisit
    # if hostname in ["stampede", "sherlock", "comet"]:
    #     fid_qsub_MPI_ONEBATCH.write( '#!/bin/bash\n' )
    #     job_name = (basename(CWD)).replace( '/', '_' )
    #     fid_qsub_MPI_ONEBATCH.write( '#SBATCH -J %s\n' % job_name )
    #     fid_qsub_MPI_ONEBATCH.write( '#SBATCH -o %s.o%%j\n' % job_name )
    #     queue = 'normal'
    #     if hostname in ['comet']:
    #         queue = 'compute'
    #     if development:
    #         queue = 'development'
    #     fid_qsub_MPI_ONEBATCH.write( '#SBATCH -p %s\n' % queue)
    #     if development:
    #         fid_qsub_MPI_ONEBATCH.write( '#SBATCH -t 00:10:00\n' )
    #     else:
    #         fid_qsub_MPI_ONEBATCH.write( '#SBATCH -t %d:00:00\n' % nhours )
    #     #fid_qsub_MPI_ONEBATCH.write( '#SBATCH --mail-user=rhiju@stanford.edu\n' )
    #     #fid_qsub_MPI_ONEBATCH.write( '#SBATCH --mail-type=ALL\n' )
    #     fid_qsub_MPI_ONEBATCH.write( '#SBATCH -n %d\n' % tot_jobs )
    #     fid_qsub_MPI_ONEBATCH.write( '#SBATCH -N %d\n' % tot_nodes )
    #     if account: fid_qsub_MPI_ONEBATCH.write( '#SBATCH -A %s\n' % account )
    #     #fid_qsub_MPI_ONEBATCH.write( 'echo $SLURM_NODELIST > nodefile.txt\n' )
    #     #fid_qsub_MPI_ONEBATCH.write( 'echo $SLURM_JOB_CPUS_PER_NODE > ncpus_per_node.txt\n' )
    #     fid_qsub_MPI_ONEBATCH.write( 'pp_jobsub.py %s -cluster_name %s -nodelist $SLURM_NODELIST -job_cpus_per_node $SLURM_JOB_CPUS_PER_NODE\n'
    #                                  % (job_file_MPI_ONEBATCH, hostname) )
    # else:
    #     fid_qsub_MPI_ONEBATCH.write( '#!/bin/bash 	 \n')
    #     fid_qsub_MPI_ONEBATCH.write( '#$ -V 	#Inherit the submission environment\n')
    #     fid_qsub_MPI_ONEBATCH.write( '#$ -cwd 	# Start job in submission directory\n')
    #     fid_qsub_MPI_ONEBATCH.write( '#$ -N %s 	# Job Name\n' % (CWD + '/' + outdir).replace( '/', '_' ) )
    #     fid_qsub_MPI_ONEBATCH.write( '#$ -j y 	# Combine stderr and stdout\n')
    #     fid_qsub_MPI_ONEBATCH.write( '#$ -o $JOB_NAME.o$JOB_ID 	# Name of the output file\n')
    #     fid_qsub_MPI_ONEBATCH.write( '#$ -pe %dway %d 	# Requests X (=12) tasks/node, Y (=12) cores total (Y must be multiples of 12, set X to 12 for lonestar)\n' % (tasks_per_node_MPI, N_MPIJOBS_ONEBATCH) )
    #     fid_qsub_MPI_ONEBATCH.write( '#$ -q normal 	# Queue name normal\n')
    if not hostname in ["stampede", "sherlock", "comet"]:
        fid_queue_MPI_ONEBATCH.write( '#!/bin/bash    \n')
        fid_queue_MPI_ONEBATCH.write( '#$ -V     #Inherit the submission environment\n')
        fid_queue_MPI_ONEBATCH.write( '#$ -cwd   # Start job in submission directory\n')
        fid_queue_MPI_ONEBATCH.write( '#$ -N %s  # Job Name\n' % (CWD + '/' + outdir).replace( '/', '_' ) )
        fid_queue_MPI_ONEBATCH.write( '#$ -j y   # Combine stderr and stdout\n')
        fid_queue_MPI_ONEBATCH.write( '#$ -o $JOB_NAME.o$JOB_ID  # Name of the output file\n')
        fid_queue_MPI_ONEBATCH.write( '#$ -pe %dway %d   # Requests X (=12) tasks/node, Y (=12) cores total (Y must be multiples of 12, set X to 12 for lonestar)\n' % (tasks_per_node_MPI, N_MPIJOBS_ONEBATCH) )
        fid_queue_MPI_ONEBATCH.write( '#$ -q normal  # Queue name normal\n')
        if nhours == 0: # for testing
            fid_queue_MPI_ONEBATCH.write( '#$ -l h_rt=00:01:00   # Run time (hh:mm:ss)\n' )
        else:
            fid_queue_MPI_ONEBATCH.write( '#$ -l h_rt=%2d:00:00  # Run time (hh:mm:ss)\n' % nhours)
        #fid_queue_MPI_ONEBATCH.write( '#$ -M rhiju@stanford.edu # Address for email notification\n')
        fid_queue_MPI_ONEBATCH.write( '#$ -m be  # Email at Begin and End of job\n')
        fid_queue_MPI_ONEBATCH.write( 'set -x    # Echo commands, use set echo with csh\n')
        fid_queue_MPI_ONEBATCH.write( 'ibrun mpi_simple_job_submit.py %s # Run the MPI python\n' % job_file_MPI_ONEBATCH)


fid.close()
fid_condor.close()
fid_qsub.close()
fid_sbatch.close()
fid_all_commands.close()

if DO_MPI:
    fid_queue_MPI.close()
    fid_queue_MPI_ONEBATCH.close()
    fid_job_MPI_ONEBATCH.close()

<<<<<<< HEAD
if len( hostname ) == 0 and bsub_file != '/dev/null':
    print 'Created bsub submission file ',bsub_file,' with ',tot_jobs, ' jobs queued. To run, type: '
    print '>source',bsub_file
    print
=======
if len( hostname ) == 0:
    print('Created bsub submission file ',bsub_file,' with ',tot_jobs, ' jobs queued. To run, type: ')
    print('>source',bsub_file)
    print()
>>>>>>> master

if hostname == 'ade':
    print('Created condor submission file ',condor_file,' with ',tot_jobs, ' jobs queued. To run, type: ')
    print('>condor_submit',condor_file)
    print()

    print('Also created bash file with all commands ',condor_file,' with ',tot_jobs, ' jobs queued. To run, type: ')
    print('>bash ', all_commands_file)
    print()

<<<<<<< HEAD
if queue_cmd == 'qsub':
    print 'Created qsub submission files ',qsub_file,' with ',tot_jobs, ' jobs queued. To run, type: '
    print '>source ',qsub_file
    print
=======
if len( hostname ) == 0:
    print('Created qsub submission files ',qsub_file,' with ',tot_jobs, ' jobs queued. To run, type: ')
    print('>source ',qsub_file)
    print()
>>>>>>> master

if queue_cmd == 'sbatch':
    print('Created sbatch submission files ',sbatch_file,' with ',tot_jobs, ' jobs queued. To run, type: ')
    print('>source ',sbatch_file)
    print()


if DO_MPI:
    if len( hostname ) == 0:
        print('Created MPI_ONEBATCH qsub submission files ',queue_file_MPI_ONEBATCH,' with ',tot_jobs, ' jobs queued. To run, type: ')
        print('>qsub ',queue_file_MPI_ONEBATCH)
        print()

    print('Created MPI submission files ',queue_file_MPI,' with ',tot_nodes, ' batches queued. To run, type: ')
    print('>source ',queue_file_MPI)

