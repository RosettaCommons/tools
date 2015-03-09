#!/usr/bin/env python

import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import sys
import math

if len( sys.argv ) < 3:
    print sys.argv[0], " <in.table> <out.png>"
    print
    print " Use to plot tables that are output by nucleobase_sample_around"
    print " application. See also Rosetta/demos/public/nucleobase_sample_around."
    exit()

#Load in the nucleobase geometry
pdb_file = open("START.pdb")

coords = []
coords_ploted = []

for line in pdb_file :
    if line[0:4] != 'ATOM' :
        continue
    atom_type = line[13]
    x = - float( line.split() [6] )
    y =   float( line.split() [7] )
    z = - float( line.split() [8] )
    coords.append( [atom_type, x, y, z] )
    if 'xz' in sys.argv[1] :
        coords_ploted.append( [x, z] )
    elif 'yz' in sys.argv[1] :
        coords_ploted.append( [y, z] )
    else :
        coords_ploted.append( [x, y] )

for i, coord1 in enumerate( coords ) :
    color = ''
    markersize = 8

    atom_type = coord1[0]
    x = coords_ploted[i][0]
    y = coords_ploted[i][1]

    if atom_type == 'H' :
        color = 'wo'
        markersize = 5
    elif atom_type == 'N'  :
        color = 'bo'
    elif atom_type == 'C'  :
        color = 'ko'
    elif atom_type == 'O'  :
        color = 'ro'
    scat = plt.plot(x, y, color)
    plt.setp(scat, 'markersize', markersize)

    for j, coord2 in enumerate( coords ) :
        if i == j :
            continue
        dist = ( (coord1[1] - coord2[1]) * (coord1[1] - coord2[1]) +
                 (coord1[2] - coord2[2]) * (coord1[2] - coord2[2]) +
                 (coord1[3] - coord2[3]) * (coord1[3] - coord2[3]) )
        dist = math.sqrt(dist)
        if dist < 1.6 :
            x_line = [x, coords_ploted[j][0]]
            y_line = [y, coords_ploted[j][1]]
            plt.plot(x_line, y_line, 'k-')


scores_file = open( sys.argv[1] )

scores = []
for line in scores_file :
    score_1d = []
    for elem in line.split() :
        score = float( elem )
        if math.isnan(score) :
            score = 1000
        score_1d.append( score )
    scores.append(score_1d)

ns = len(scores)
scores = np.array(scores)

sum_bound = 0
n_data = 0
for i in xrange(ns):
    for j in xrange(ns):
        if i == 0 or i == ns - 1 or j == 0 or j == ns - 1 :
            sum_bound += scores[i,j]
            n_data += 1

baseline = sum_bound / n_data

for i in xrange(ns):
    for j in xrange(ns):
        if math.isnan(scores[i,j]) :
            print i, ' ', j
        scores[i,j] -= baseline

min_score = scores.min()

sample_grid = 0.2
sample_max = (len( scores ) - 1) * 0.5 * sample_grid
X = np.arange(-sample_max, sample_max + sample_grid, sample_grid)
Y = X

grid_size = abs(min_score) / 100.0
levels = np.arange(min_score, grid_size*40, grid_size)

CS = plt.contourf(X, Y, scores, levels)

#plt.clabel(CS, inline=1, fontsize=10)
plt.axis([-10,10,-10,10])
plt.colorbar()


if len(sys.argv) == 3:
    plt.savefig( sys.argv[2], dpi=200 )
else :
    plt.show()

