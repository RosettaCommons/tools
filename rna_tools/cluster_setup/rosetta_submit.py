#!/usr/bin/python

from sys import argv,exit
from os import system,getcwd,popen
from os.path import basename,dirname,expanduser,exists
import string

def Help():
    print argv[0]+' <text file with rosetta command> <outdir> <# jobs>  [# hours]'
    exit()

if len( argv ) < 4:
    Help()

infile = argv[1]
outdir = argv[2]
try:
    n_jobs = int( argv[3] )
except:
    print 'NEED TO SUPPLY NUMBER OF JOBS'


tasks_per_node_MPI = 12 # lonestar

hostname = ''
hostname_tag = popen( 'hostname' ).readlines()[0]
if hostname_tag.find( 'ls4' ) > -1: hostname = 'lonestar'
if hostname_tag.find( 'DasLab' ) > -1: hostname = 'ade'

save_logs = False
if argv.count( '-save_logs' )>0:
    save_logs = True
    pos = argv.index( '-save_logs' )
    del( argv[pos] )

nhours = 16
if len( argv ) > 4:
    nhours = int( argv[4] )
    if ( nhours > 168 ):  Help()

lines = open(infile).readlines()

bsub_file = 'bsubMINI'
condor_file = 'condorMINI'
qsub_file = 'qsubMINI'
qsub_file_MPI = 'qsubMPI'
qsub_file_MPI_ONEBATCH = 'qsubMPI_ONEBATCH'
job_file_MPI_ONEBATCH = 'MPI_ONEBATCH.job'

fid = open( bsub_file,'w')
fid_condor = open( condor_file,'w')
fid_qsub = open( qsub_file,'w')

fid_qsub_MPI = open( qsub_file_MPI,'w')
fid_job_MPI_ONEBATCH = open( job_file_MPI_ONEBATCH,'w')
fid_qsub_MPI_ONEBATCH = open( qsub_file_MPI_ONEBATCH,'w')

tot_jobs = 0

universe = 'vanilla';
fid_condor.write('+TGProject = "TG-MCB090153"\n')
fid_condor.write('universe = %s\n' % universe)
fid_condor.write('notification = never\n')

HOMEDIR = expanduser('~')
CWD = getcwd()

qsub_file_dir = 'qsub_files/'
if not exists( qsub_file_dir ): system( 'mkdir '+qsub_file_dir )

qsub_file_dir_MPI = 'qsub_files_MPI/'
if not exists( qsub_file_dir_MPI ): system( 'mkdir '+qsub_file_dir_MPI )

command_lines_explicit = []

for line in  lines:

    if len(line) == 0: continue
    if line[0] == '#': continue
    #if string.split( line[0]) == []: continue

    dir = outdir + '/$(Process)/'
    command_line = line[:-1].replace( 'out:file:silent  ','out:file:silent ').replace( '-out:file:silent ', '-out:file:silent '+dir)
    command_line = command_line.replace( '-out::file::silent ', '-out::file::silent '+dir)
    command_line = command_line.replace( '-out:file:o ', '-out:file:o '+dir)
    command_line = command_line.replace( '-o ', '-o '+dir)
    command_line = command_line.replace( '-seed_offset 0', '-seed_offset $(Process)')
    command_line = command_line.replace( 'macosgcc', 'linuxgcc')
    command_line = command_line.replace( 'Users', 'home')
    command_line = command_line.replace( '~/', HOMEDIR+'/')
    command_line = command_line.replace( '/home/rhiju',HOMEDIR)

    cols = string.split( command_line )

    if len( cols ) == 0: continue

    if '-total_jobs' in cols:
        pos = cols.index( '-total_jobs' )
        cols[ pos+1 ] = '%d' % n_jobs
        command_line = string.join( cols )
    if '-job_number' in cols:
        pos = cols.index( '-job_number' )
        cols[ pos+1 ] = '$(Process)'
        command_line = string.join( cols )

    if save_logs:
        outfile_general = '$(Process).out'
        errfile_general = '$(Process).err'
    else:
        outfile_general = '/dev/null'
        errfile_general = '/dev/null'

    for i in range( n_jobs ):
        dir_actual = dir.replace( '$(Process)', '%d' % i)
        system( 'mkdir -p '+ dirname(dir_actual) )

        outfile = outfile_general.replace( '$(Process)', '%d' % i )
        errfile = errfile_general.replace( '$(Process)', '%d' % i )

        command =  'bsub -W %d:0 -o %s -e %s ' % (nhours, outfile, errfile )
        command_line_explicit = command_line.replace( '$(Process)', '%d' % i )
        command += command_line_explicit
        fid.write( command + '\n')

        # qsub
        qsub_submit_file = '%s/qsub%d.sh' % (qsub_file_dir, tot_jobs )
        fid_qsub_submit_file = open( qsub_submit_file, 'w' )
        fid_qsub_submit_file.write( '#!/bin/bash\n'  )
        fid_qsub_submit_file.write('#PBS -N %s\n' %  (CWD+'/'+dir_actual[:-1]).replace( '/', '_' ) )
        fid_qsub_submit_file.write('#PBS -o %s\n' % outfile)
        fid_qsub_submit_file.write('#PBS -e %s\n' % errfile)
        fid_qsub_submit_file.write('#PBS -l walltime=48:00:00\n\n')
        fid_qsub_submit_file.write( 'cd %s\n\n' % CWD )
        fid_qsub_submit_file.write( command_line_explicit+' > /dev/null 2> /dev/null \n' )
        fid_qsub_submit_file.close()

        fid_qsub.write( 'qsub %s\n' % qsub_submit_file )


        # MPI job file
        fid_job_MPI_ONEBATCH.write( '%s ;;; %s\n' % (CWD, command_line_explicit) )

        command_lines_explicit.append( command_line_explicit )
        tot_jobs += 1

    EXE = cols[ 0 ]
    if not exists( EXE ): EXE = EXE.replace( 'linux', 'macos' )
    if not exists( EXE ): EXE = EXE.replace( 'macos', 'linux' )
    if not exists( EXE ):
        EXE = HOMEDIR + '/src/rosetta_TRUNK/rosetta_source/bin/'+EXE
    if not exists( EXE ):
        EXE = HOMEDIR + '/src/mini/bin/'+EXE
        assert( exists( EXE ) )
    arguments = string.join( cols[ 1: ] )


    fid_condor.write('\nexecutable = %s\n' % EXE )
    fid_condor.write('arguments = %s\n' % arguments)
    if save_logs:
        fid_condor.write( 'output = %s\n' % outfile_general )
        fid_condor.write( 'error  = %s\n' % errfile_general )
    fid_condor.write('Queue %d\n' % n_jobs )


N_MPIJOBS_ONEBATCH =  tot_jobs
tot_nodes = ( tot_jobs / tasks_per_node_MPI )

if ( N_MPIJOBS_ONEBATCH % tasks_per_node_MPI != 0 ):  # checking modulo 12 (or whatever the number of cores/node)
    N_MPIJOBS_ONEBATCH += ( tot_jobs/ tasks_per_node_MPI + 1) * tasks_per_node_MPI
    tot_nodes += 1


count = 0
for n in range( tot_nodes ):

    job_submit_file_MPI = '%s/qsubMPI%d.job' % (qsub_file_dir_MPI, n )
    fid_job_submit_file_MPI  = open( job_submit_file_MPI, 'w' )

    for m in range( tasks_per_node_MPI ):
        count = count + 1
        if ( count <= tot_jobs ):
            command_line_explicit = command_lines_explicit[ count-1 ]
            fid_job_submit_file_MPI.write( '%s ;;; %s\n' % (CWD, command_line_explicit) )

    fid_job_submit_file_MPI.close()

    # qsub MPI
    jobname= (CWD + '/' + outdir).replace( '/', '_' )
    jobname = jobname[-30:]
    qsub_submit_file_MPI = '%s/qsubMPI%d.sh' % (qsub_file_dir_MPI, n )
    fid_qsub_submit_file_MPI = open( qsub_submit_file_MPI, 'w' )
    fid_qsub_submit_file_MPI.write( '#!/bin/bash 	 \n')
    fid_qsub_submit_file_MPI.write( '#$ -V 	#Inherit the submission environment\n')
    fid_qsub_submit_file_MPI.write( '#$ -cwd 	# Start job in submission directory\n')
    fid_qsub_submit_file_MPI.write( '#$ -N %s 	# Job Name\n' % jobname )
    fid_qsub_submit_file_MPI.write( '#$ -j y 	# Combine stderr and stdout\n')
    fid_qsub_submit_file_MPI.write( '#$ -o $JOB_NAME.o$JOB_ID 	# Name of the output file\n')
    fid_qsub_submit_file_MPI.write( '#$ -pe %dway %d 	# Requests X (=12) tasks/node, Y (=12) cores total (Y must be multiples of 12, set X to 12 for lonestar)\n' % (tasks_per_node_MPI, tasks_per_node_MPI) )
    fid_qsub_submit_file_MPI.write( '#$ -q normal 	# Queue name normal\n')
    if nhours == 0: # for testing
        fid_qsub_submit_file_MPI.write( '#$ -l h_rt=00:01:00 	# Run time (hh:mm:ss)\n' )
    else:
        fid_qsub_submit_file_MPI.write( '#$ -l h_rt=%02d:00:00 	# Run time (hh:mm:ss)\n' % nhours)
    #fid_qsub_submit_file_MPI.write( '#$ -M rhiju@stanford.edu	# Address for email notification\n')
    fid_qsub_submit_file_MPI.write( '#$ -m be 	# Email at Begin and End of job\n')
    fid_qsub_submit_file_MPI.write( 'set -x 	# Echo commands, use set echo with csh\n')
    fid_qsub_submit_file_MPI.write( 'ibrun mpi_simple_job_submit.py %s	# Run the MPI python\n' % job_submit_file_MPI)
    fid_qsub_submit_file_MPI.close()


    fid_qsub_MPI.write( 'qsub %s\n' % qsub_submit_file_MPI )



fid_qsub_MPI_ONEBATCH.write( '#!/bin/bash 	 \n')
fid_qsub_MPI_ONEBATCH.write( '#$ -V 	#Inherit the submission environment\n')
fid_qsub_MPI_ONEBATCH.write( '#$ -cwd 	# Start job in submission directory\n')
fid_qsub_MPI_ONEBATCH.write( '#$ -N %s 	# Job Name\n' % (CWD + '/' + outdir).replace( '/', '_' ) )
fid_qsub_MPI_ONEBATCH.write( '#$ -j y 	# Combine stderr and stdout\n')
fid_qsub_MPI_ONEBATCH.write( '#$ -o $JOB_NAME.o$JOB_ID 	# Name of the output file\n')
fid_qsub_MPI_ONEBATCH.write( '#$ -pe %dway %d 	# Requests X (=12) tasks/node, Y (=12) cores total (Y must be multiples of 12, set X to 12 for lonestar)\n' % (tasks_per_node_MPI, N_MPIJOBS_ONEBATCH) )
fid_qsub_MPI_ONEBATCH.write( '#$ -q normal 	# Queue name normal\n')
if nhours == 0: # for testing
    fid_qsub_MPI_ONEBATCH.write( '#$ -l h_rt=00:01:00 	# Run time (hh:mm:ss)\n' )
else:
    fid_qsub_MPI_ONEBATCH.write( '#$ -l h_rt=%2d:00:00 	# Run time (hh:mm:ss)\n' % nhours)
#fid_qsub_MPI_ONEBATCH.write( '#$ -M rhiju@stanford.edu	# Address for email notification\n')
fid_qsub_MPI_ONEBATCH.write( '#$ -m be 	# Email at Begin and End of job\n')
fid_qsub_MPI_ONEBATCH.write( 'set -x 	# Echo commands, use set echo with csh\n')
fid_qsub_MPI_ONEBATCH.write( 'ibrun mpi_simple_job_submit.py %s	# Run the MPI python\n' % job_file_MPI_ONEBATCH)


fid.close()
fid_condor.close()
fid_qsub.close()
fid_qsub_MPI.close()
fid_qsub_MPI_ONEBATCH.close()
fid_job_MPI_ONEBATCH.close()


if len( hostname ) == 0:
    print 'Created bsub submission file ',bsub_file,' with ',tot_jobs, ' jobs queued. To run, type: '
    print '>source',bsub_file
    print

if len( hostname ) == 0 or hostname == 'ade':
    print 'Created condor submission file ',condor_file,' with ',tot_jobs, ' jobs queued. To run, type: '
    print '>condor_submit',condor_file
    print

if len( hostname ) == 0:
    print 'Created qsub submission files ',qsub_file,' with ',tot_jobs, ' jobs queued. To run, type: '
    print '>source ',qsub_file
    print

if len( hostname ) == 0:
    print 'Created MPI_ONEBATCH qsub submission files ',qsub_file_MPI_ONEBATCH,' with ',tot_jobs, ' jobs queued. To run, type: '
    print '>qsub ',qsub_file_MPI_ONEBATCH
    print

if len( hostname ) == 0 or hostname == 'lonestar':
    print 'Created MPI submission files ',qsub_file_MPI,' with ',tot_jobs, ' jobs queued. To run, type: '
    print '>source ',qsub_file_MPI
