#!/usr/bin/env python2.7-mpi

from mpi4py import MPI

from sys import argv,exit,stdout
import string
from os.path import basename,exists
from time import sleep

######################################################################
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options, option_name_exist, get_option_name_args_safe

from DAG_continuous import *
######################################################################

#####!/usr/bin/env python2.7-mpi


#mpiexec -n 4 python2.7-mpi $WORK/SWA_RNA_python/SWA_dagman_python/dagman/DAG_continuous_MPI_wrapper.py -num_slave_jobs 3 -dagman_file rna_build.dag  

#mpiexec -n 4 python $WORK/SWA_RNA_python/SWA_dagman_python/dagman/DAG_continuous_MPI_wrapper.py -num_slave_jobs 3 -dagman_file rna_build.dag  

#mpiexec -n 4 ~/SWA_RNA_python/SWA_dagman_python/dagman/DAG_continuous_MPI_wrapper.py -num_slave_jobs 3 -dagman_file rna_build.dag  

#aprun -n 4  python2.7-mpi /lustre/scratch/fcchou/SWA_RNA_python/SWA_dagman_python/dagman/DAG_continuous_MPI_wrapper.py -num_slave_jobs 3 -dagman_file rna_build.dag  

#mpiexec -n 4 /lustre/scratch/fcchou/SWA_RNA_python/SWA_dagman_python/dagman/DAG_continuous_MPI_wrapper.py -num_slave_jobs 3 -dagman_file rna_build.dag  

#aprun -n 4  python2.7-mpi  /lustre/scratch/fcchou/mpi4py_lib/test_script/helloworld.py 

#ibrun  /work/01462/parins/SWA_RNA_python/SWA_dagman_python/dagman/DAG_continuous_MPI_wrapper.py -num_slave_jobs 23 -dagman_file rna_build.dag >MPI_log.out 2>MPI_log.err

#> MPI_LOG.out 2> MPI_LOG.err

######################################################################

MPI_comm = MPI.COMM_WORLD
MPI_rank = MPI_comm.Get_rank() #Numbering of current core
MPI_size = MPI_comm.Get_size() #Total number of cores


START_argv=copy.deepcopy(argv)

DAG_continuous_args=list_to_string(START_argv[1:])

if(option_name_exist(DAG_continuous_args, "slave_wall_time")):
	error_exit_with_message("slave_wall_time option passed into DAG_continous_MPI_wrapper.py!")

if(option_name_exist(DAG_continuous_args, "num_slave_jobs")==False):
	error_exit_with_message("num_slave_jobs option was not passed into DAG_continous_MPI_wrapper.py!")	

num_slave_jobs=get_option_name_args_safe(DAG_continuous_args, "num_slave_jobs", 0)

#Consistency check.
if(num_slave_jobs!=(MPI_size-1)): error_exit_with_message("num_slave_jobs=(%s)!=(%s)=(MPI_size-1)" %(num_slave_jobs,MPI_size-1))



if(MPI_rank == 0): #master_job

	print "Enter DAG_continuous_MPI_wrapper.py: master process"

	out_log="master_log.out"
	err_log="master_log.err"

	if(exists(out_log)): submit_subprocess("rm %s" %(out_log))
	if(exists(err_log)): submit_subprocess("rm %s" %(err_log))

	errfile= safe_open(err_log, "w")
	outfile= safe_open(out_log, "w")

	os.dup2(outfile.fileno(), sys.stdout.fileno())
	os.dup2(errfile.fileno(), sys.stderr.fileno())

	DAG_continuous( ("DAG_continous.py %s " %(DAG_continuous_args)).split() )

	print "Exit DAG_continuous_MPI_wrapper.py: master process"

else: #slave_jobs

	print "Enter DAG_continuous_MPI_wrapper.py: slave process %s of %s" % (MPI_rank, MPI_size-1)

	#slave_wall_time=parse_options( argv, "slave_wall_time", 48) #In MPI mode, slave_wall_time is set when submitting the MPI job.

	slave_job_ID=MPI_rank-1 #Start from 0 to num_slave_jobs-1 (or equivalently MPI_size-2)

	MPI_DAG_slave_wrapper(slave_job_ID, num_slave_jobs, SLAVE_EXE)


	print "Exit DAG_continuous_MPI_wrapper.py: slave process %s of %s" % (MPI_rank, MPI_size-1)

