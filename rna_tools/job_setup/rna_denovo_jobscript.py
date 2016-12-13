#!/usr/bin/env python
import argparse
import os
import os.path as path
import shutil


parser = argparse.ArgumentParser(
    description='Setup rna_denovo jobscripts for parallel runs.')
parser.add_argument(
    'job_file', help='Roseta command lines')
parser.add_argument(
    '-proc', help='Number of processes (default 100)', default=100, type=int)
parser.add_argument(
    '-out_folder', help='Folder for the output files', default='out')
parser.add_argument(
    '-out_script', help='Job Script for submission',
    default='rna_denovo.joblist')

args = parser.parse_args()

# Get and parse rna_denovo cmdline
with open(args.job_file) as script_f:
    orig_cmdline = script_f.read().strip()
split_cmd = orig_cmdline.split()
for i, elem in enumerate(split_cmd):
    if elem == '-out:file:silent':
        out_file_silent = split_cmd[i + 1]
        new_cmd = ' '.join(split_cmd[:i] + split_cmd[(i + 2):])
        break
else:
    new_cmd = orig_cmdline
    out_file_silent = 'default.out'


if path.isdir(args.out_folder):
    shutil.rmtree(args.out_folder)
os.mkdir(args.out_folder)

out_script = open(args.out_script, 'w')
for i in xrange(args.proc):
    folder = path.join(args.out_folder, str(i))
    os.mkdir(folder)
    out_file = path.join(folder, out_file_silent)
    cmdline = new_cmd + ' -out:file:silent %s' % out_file
    out_script.write('./\t%s\n' % cmdline)
out_script.close()
