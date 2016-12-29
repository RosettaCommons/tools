#!/usr/bin/env python
import argparse
from parse_options import get_ints

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

parser = argparse.ArgumentParser(description='Sort a slient file according to its scores and output selected ones')
parser.add_argument('file', help='Input silent file')
parser.add_argument('-select', nargs='+', help='Selected output file index ranked by scores, ex. 1-9 for lowest score 9 decoys')
parser.add_argument('-o', default='silent_sort.out', help='Filename of output silent file')
parser.add_argument('-term', default='score', help='Option to input a specific column to sort by (e.g. fa_atr, rms); default is score')
args = parser.parse_args()

select_idx = None
if args.select is not None:
    select_idx = []
    for i in args.select:
        get_ints(i, select_idx)

#Reads the silent file
header = ''
scores_data = []
new_score_data = None
is_header_over = False
for line in open(args.file):
    if not is_header_over:
        if len(line) > 5 and line[:5] == 'SCORE':
            if is_number(line.split()[1]):
                is_header_over = True
            else:
                term_col = line.split().index(args.term)
    if not is_header_over:
        header += line
    else:
        if len(line) > 5 and line[:5] == 'SCORE':
            if new_score_data is not None:
                scores_data.append(new_score_data)
            new_score_data = [0, '']
            term = float(line.split()[term_col])
            new_score_data[0] = term
            new_score_data[1] = line
        else:
            new_score_data[1] += line
else:
    scores_data.append(new_score_data)
#Sort the data
sorted_data = sorted(scores_data, key = lambda x : x[0])

#Output
out = open(args.o, 'w')
out.write(header)
if select_idx is None: #Output everything
    for data in sorted_data:
        out.write(data[1])
else: #Output selected data
    for i in select_idx:
        out.write(sorted_data[i-1][1])
out.close()
