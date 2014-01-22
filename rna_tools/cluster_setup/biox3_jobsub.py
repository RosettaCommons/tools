#!/usr/bin/env python
import qsub_torque
import sys
import os
import os.path

doc_string = '''
Usage: biox3_jobsub.py <job_script>
'''
print doc_string


def load_jobfile(filename):
    jobs = []
    for line in open(filename):
        work_dir, cmdline = line.strip().split('\t', 1)
        if work_dir != '' and cmdline != '':
            jobs.append((os.path.abspath(work_dir), cmdline))
    return jobs

cwd = os.getcwd()

for path, cmdline in load_jobfile(sys.argv[1]):
    os.chdir(path)
    for i in xrange(50):
        try:
            qsub_torque.run(cmdline)
            break
        except:
            continue
    else:
        raise Exception("Fail to submit!")
    os.chdir(cwd)
