import random
import subprocess
import pp
from os.path import abspath
import os
import time


def parse_nodefile(inp):
    '''
    Parse nodelist from file. Currently requires SLURM_NODELIST format.
    '''
    inp.replace('\n', ',')
    node_list = []
    while len(inp) != 0:
        if inp[0] == ',':
            inp = inp[1:]
        nodes = inp.split(', ', 1)
        if '[' not in nodes[0]:
            node_list.append(nodes[0])
            if len(nodes) > 1:
                inp = nodes[1]
            else:
                inp = ''
        else:
            start_brc = inp.index('[')
            end_brc = inp.index(']')
            initial = inp[:start_brc]
            nodeids = inp[(start_brc+1):end_brc]
            for i in nodeids.split(','):
                if '-' not in i:
                    node_list.append(initial + i)
                else:
                    start_id = int(i.split('-')[0])
                    end_id = int(i.split('-')[1])
                    for j in xrange(start_id, end_id+1):
                        node_list.append(initial + '%03d' % j)
            inp = inp[(end_brc+1):]
    return node_list


def jobserver_init(cluster_name, nodelist = '', job_cpus_per_node=''):
    if cluster_name in ['stampede']:
        submit_command = 'ibrun'
    else:
        submit_command = 'srun'
    port = 32568
    socket_timeout = 3600 * 24
    key_phrase = '%x' % random.randrange(256**5)

    if len( nodelist ) == 0:
        nodelist = open( 'nodefile.txt' ).read().strip()
    if len( job_cpus_per_node ) == 0:
        job_cpus_per_node = open('ncpus_per_node.txt').read().strip()
    nodes = parse_nodefile(nodelist)
    nworkers = int(job_cpus_per_node.split('(')[0])
    ncpus = nworkers * len(nodes)

    print( 'Submitting the ppserver...' )
    submit_cmdline = [submit_command, '/home/rhiju/src/pp-1.6.4/ppserver.py']
    if not nworkers is None:
        submit_cmdline += ['-w', str(nworkers)]
    if not socket_timeout is None:
        submit_cmdline += ['-k', str(socket_timeout)]
    if not port is None:
        submit_cmdline += ['-p', str(port)]
    if not key_phrase is None:
        submit_cmdline += ['-s', key_phrase]
    subprocess.Popen(submit_cmdline)
    if port is None:
        ppservers = tuple([node for node in nodes])
    else:
        ppservers = tuple([node + ':' + str(port) for node in nodes])
    jobserver = pp.Server(
        ncpus=ncpus, ppservers=ppservers, secret=key_phrase,
        socket_timeout=socket_timeout)
    time.sleep(30)
    return jobserver, ncpus


def load_jobfile(filename):
    work_dir_list = []
    cmdline_list = []
    for line in open(filename):
        try:
            work_dir, cmdline = line.strip().split('\t', 1)
            if work_dir != '' and cmdline != '':
                work_dir_list.append(abspath(work_dir))
                cmdline_list.append(cmdline)
        except:
            continue
    print( cmdline_list )
    return work_dir_list, cmdline_list


def submit_cmdline(work_dir, cmdline):
    cwd = os.getcwd()
    os.chdir(work_dir)
    output = ''
    returncode = 0
    try:
        output, error = subprocess.Popen(
            cmdline.split(),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT ).communicate()
        if error:
            print( error )
            output = error
            returncode = 1
    except subprocess.CalledProcessError as err:
        output = err.output
        returncode = err.returncode
    os.chdir(cwd)
    return output, returncode
