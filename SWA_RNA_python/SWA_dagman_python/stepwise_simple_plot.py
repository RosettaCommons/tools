#!/usr/bin/python

from matplotlib import pyplot as plt 


silent_file='region_FINAL.out'
xvar = 'all_rms'
yvar = 'score'

score_keys = [line for line in open(silent_file, 'r').readlines() if 'SCORE' in line and 'description' in line][0].split()
score_lines = [line for line in open(silent_file, 'r').readlines() if 'SCORE' in line and 'description' not in line]

xdata = [float(cols[score_keys.index(xvar)]) for cols in [line.split() for line in score_lines]]
ydata = [float(cols[score_keys.index(yvar)]) for cols in [line.split() for line in score_lines]]

#for s,r in zip(score_,rms_):
#	print s,', ',r 

plot = plt.plot(xdata,ydata, marker='o', color='r', linestyle=' ')
plt.xlabel(xvar)
plt.ylabel(yvar)
plt.title(xvar+' vs. '+yvar)
plt.show()
