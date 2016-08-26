#!/usr/bin/env python
import argparse
import os
import os.path as path
import shutil
from warnings import warn
from rosetta_exe import rosetta_exe


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

parser = argparse.ArgumentParser(
    description='Setup parallel minimize for a silent file'
)
parser.add_argument('-silent', help='Input silent file', required=True)
parser.add_argument(
    '-tag', help='Outfile tag', default='MIN'
)
parser.add_argument(
    '-proc', help='Number of processes (default 0)', default=20, type=int
)
parser.add_argument(
    '-nstruct', help='Number of structures to be minimzie (default 1000)',
    default=1000, type=int
)
parser.add_argument(
    '-out_folder', help='Folder for the output files',
    default='out'
)
parser.add_argument(
    '-out_script', help='Output script file for cmdlines',
    default='min_cmdlines'
)
parser.add_argument(
    'extra_opt', help='Extra options for Rosetta'
)

args = parser.parse_args()
if args.proc > args.nstruct:
    raise Exception("More proc than nstruct! Wasting resource...")
if not path.isfile(args.silent):
    raise Exception("Cannot locate the input silent file!")


#Reads the silent file
header = ''
scores_data = []
new_score_data = None
is_header_over = False
for line in open(args.silent):
    if not is_header_over:
        if (len(line) > 5 and line[:5] == 'SCORE'
                and is_number(line.split()[1])):
            is_header_over = True
    if not is_header_over:
        header += line
    else:
        if len(line) > 5 and line[:5] == 'SCORE':
            if new_score_data is not None:
                scores_data.append(new_score_data)
            new_score_data = [0, bytearray()]
            score = float(line.split()[1])
            new_score_data[0] = score
            new_score_data[1] = line
        else:
            new_score_data[1] += line
else:
    scores_data.append(new_score_data)
#Sort the data
#sorted_data = sorted(scores_data, key=lambda x: x[0])
scores_data.sort(key=lambda x: x[0])
n_data = len(scores_data)
if n_data < args.nstruct:
    warn("Warning: n_data in silent file smaller than nstruct!")
    nstruct = n_data
else:
    nstruct = args.nstruct

if path.isdir(args.out_folder):
    shutil.rmtree(args.out_folder)
os.mkdir(args.out_folder)

out_script = open(args.out_script, 'w')
silent_fh = []
cmdline_base = rosetta_exe('rna_minimize')
cmdline_base += ' ' + args.extra_opt + ' -in:file:silent '
for i in xrange(args.proc):
    folder = path.join(args.out_folder, str(i))
    os.mkdir(folder)
    silent_file = path.join(folder, '%d.silent' % i)
    silent_fh.append(open(silent_file, 'w'))
    silent_fh[-1].write(header)
    out_silent = path.join(folder, args.tag + '.out')
    cmdline = cmdline_base + silent_file + ' -out:file:silent %s' % out_silent
    #out_script.write('./\t%s\n' % cmdline)
    out_script.write('%s\n' % cmdline)
out_script.close()


j = 0
for i, data in enumerate(scores_data):
    if i >= nstruct:
        break
    silent_fh[j].write(data[1])
    j += 1
    if j == args.proc:
        j = 0

for i in silent_fh:
    i.close()
