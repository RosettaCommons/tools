#!/usr/bin/env python
import sys
import argparse
import pp_util


parser = argparse.ArgumentParser(
    description='Run jobs with parallel python.')
parser.add_argument('job_script', type=str, help='Job script to be submitted.')
parser.add_argument(
    '-cluster_name', type=str, help='Name of the cluster',
    metavar='str', default='stampede')
parser.add_argument(
    '-nodelist', type=str, help='node list, should provided by SLURM',
    metavar='str', default='')
parser.add_argument(
    '-job_cpus_per_node', type=str, help='job cpus per node, should provided by SLURM',
    metavar='str', default='')
args = parser.parse_args()
work_dir_list, cmdline_list = pp_util.load_jobfile(args.job_script)
if args.cluster_name in ['stampede', 'sherlock', 'comet']:
    jobserver, ncpus = pp_util.jobserver_init( args.cluster_name, args.nodelist, args.job_cpus_per_node )
else:
    raise argparse.ArgumentError("Invalid cluster_name!")
active_nodes = jobserver.get_active_nodes()
if 'local' in active_nodes:
    active_nodes.pop('local')
active_cpus = sum(active_nodes.values())
print( 'Active_nodes:', active_nodes )
print( 'N_cpus = %d, Active_cpus = %d' % (ncpus, active_cpus) )
assert (ncpus == active_cpus)

print( 'Starting submitting jobs...' )
print( '# of jobs = %d, # of cpus = %d' % (len(cmdline_list), ncpus) )
if len(cmdline_list) > ncpus:
    sys.stderr.write('WARNING: More jobs than the total amount of CPUs\n')

jobs = []
for work_dir, cmdline in zip(work_dir_list, cmdline_list):
    jobs.append(
        jobserver.submit(
            pp_util.submit_cmdline, (work_dir, cmdline),
            modules=('subprocess', 'os')))
jobserver.wait()

for i, job in enumerate(jobs):
    print( '#####Job %4d#####' % i )
    output, returncode = job()
    if returncode != 0:
        print(
            'ERROR: Job %4d returned non-zero exit status %d!!!'
            % (i, returncode))
    print( output )
print( '####################' )

jobserver.print_stats()
jobserver.destroy()
