#!/usr/bin/env python2.7-mpi
from mpi4py import MPI
import sys
from os.path import exists
import os
import subprocess

doc_string = '''
Usage: mpi4py_jobsub.py <job_script>
This code depends on a working installation of mpi4py.
Check if python2.7-mpi executable exists in your computer.
'''
print doc_string

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

print "Entering process %s of %s..." % (rank, size)

if rank == 0:
    job_file = sys.argv[1]
    assert(exists(job_file))

    working_dir_0 = ''
    cmdline_0 = ''
    for line_index, line in enumerate(open(job_file)):
        if line_index >= size:
            sys.stderr.write(
                'Error: number of jobs in job_file '
                '> number of mpi processes!\n')
            comm.Abort(1)
            assert(False)
        else:
            working_dir, cmdline = line.strip().split('\t', 1)
            if line_index == 0:
                working_dir_0 = working_dir
                cmdline_0 = cmdline
            else:
                comm.send([working_dir, cmdline], dest=line_index, tag=1)

    if line_index + 1 < size:
        sys.stderr.write(
            'Warning: %s processes initiated but only %s running jobs!!!\n'
            % (size, line_index + 1))
        for i in xrange(line_index, size):
            comm.send([], dest=i, tag=1)

    os.chdir(working_dir_0)
    err_code = subprocess.call(cmdline_0, shell=True)
    if err_code != 0:
        print(
            "Warning: Process %s failed with error code %s!!!"
            % (rank, err_code))

else:
    cmd_data = comm.recv(source=0, tag=1)
    if cmd_data != []:
        working_dir = cmd_data[0]
        cmdline = cmd_data[1]
        os.chdir(working_dir)
        err_code = subprocess.call(cmdline, shell=True)
        if err_code != 0:
            print(
                "Warning: Process %s failed with error code %s !!!"
                % (rank, err_code))

comm.Barrier()
