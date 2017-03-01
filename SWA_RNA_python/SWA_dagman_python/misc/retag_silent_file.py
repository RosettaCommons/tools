#!/usr/bin/python

#######################################################################

from os.path import exists
from sys import argv, exit
import operator
import string

#######################################################################

fname = '__empty__.txt'
score_type = 'score'

if (len(argv) > 1):  fname = argv[1]
if (len(argv) > 2):  score_type = argv[2]

while not exists(fname):
    fname = raw_input("Please enter the name of a valid silent file: ")
    if fname == 'exit':  exit(1)

#######################################################################

scores = {}
silent_structs = {}
score_idx = -1
ssid = -1
header = ''

### read in and store scores and silent structures
for line in open(fname, 'r').readlines():
    if ('SCORE' in line):
        if 'description' in line:
            while score_type not in line:
                print line
                score_type = raw_input("Please choose a score_type for reordering: ")
                if score_type == 'exit':  exit(1)
            score_idx = line.split().index(score_type)
            header += line
            continue
        else:
            assert(score_idx > -1)
            score = float(line.split()[score_idx])
            ssid += 1
            scores[ ssid ] = score
            silent_structs[ ssid ] = [ ]
    if (ssid == -1):
        header += line
    else:
        silent_structs[ ssid ].append(line.replace('empty_tag', 'S_%d' % ssid))

### sort by score_type
sorted_scores = sorted(scores.items(), key=operator.itemgetter(1))

### rewrite silent file
fout = open(fname, 'w')
fout.write(header)
for (ssid, score) in sorted_scores:
    print '%s: %f\tssid: %d' % (score_type, score, ssid)
    fout.write(string.join(silent_structs[ ssid ], ''))
    #for line in silent_structs[ ssid ]:
    #    fout.write(line)
print "Found %d silent structures." % (len(sorted_scores))

