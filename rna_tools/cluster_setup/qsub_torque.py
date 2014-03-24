#!/usr/bin/env python
import subprocess
from os import remove
import argparse


def run(
    cmd, name='default', time='36:00:00',
    ppn=1, qsub_script='job_script.qsub',
    save_script=False, output_only=False
):
    job = open(qsub_script, 'w')
    job.write('#!/bin/bash\n')
    job.write('#PBS -N %s\n' % name)
    job.write('#PBS -o $PBS_JOBNAME.out\n')
    job.write('#PBS -j oe\n')
    job.write('#PBS -l walltime=%s\n' % time)
    job.write('#PBS -l nodes=1:ppn=%d\n' % ppn)
    job.write('#PBS -V\n\n')
    job.write('cd $PBS_O_WORKDIR\n')
    job.write('%s\n' % cmd)
    job.close()

    if output_only:
        return
    else:
        subprocess.check_call(['qsub', qsub_script])
    if not save_script:
        remove(qsub_script)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Perform qsub of a command under TORQUE.')
    parser.add_argument('cmd', type=str, nargs='+', help='Command to be qsub')
    parser.add_argument(
        '-name', type=str, help='Name of the job',
        metavar='str', default='default')
    parser.add_argument(
        '-time', type=str, help="Max walltime ('13:00:00' for 13 hr)",
        metavar='str', default='36:00:00')
    parser.add_argument(
        '-ppn', type=int, help='N processor used in the job',
        metavar='int', default=1)
    args = parser.parse_args()
    run(' '.join(args.cmd), args.name, args.time, args.ppn)
