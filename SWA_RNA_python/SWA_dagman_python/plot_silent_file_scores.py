#!/usr/bin/python

from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse

parser = argparse.ArgumentParser(description='Plot scores terms from silent files.')
parser.add_argument('silent_file', help='Name of the Silent File')
parser.add_argument('out_file', help='Name of Output File', default='plot')
parser.add_argument('x', help='X-axis variable', default='rms')
parser.add_argument('y', help='Y-axis variable', default='score')
args=parser.parse_args()

silent_file=args.silent_file
out_file=args.out_file+'.pdf'
xvar=args.x
yvar=args.y
xlab=ur'RMSD (\u00c5)'#xvar.replace('_',' ').upper()
ylab='Rosetta Energy'#yvar.replace('_',' ').upper()

score_keys = [line for line in open(silent_file, 'r').readlines() if 'SCORE' in line and 'description' in line][0].split()
score_lines = [line for line in open(silent_file, 'r').readlines() if 'SCORE' in line and 'description' not in line]

xdata = [float(cols[score_keys.index(xvar)]) for cols in [line.split() for line in score_lines]]
ydata = [float(cols[score_keys.index(yvar)]) for cols in [line.split() for line in score_lines]]

#for s,r in zip(score_,rms_):
#	print s,', ',r

pp=PdfPages(out_file)
plot = plt.plot(xdata,ydata, marker='o', color='r', linestyle=' ')
plt.xlabel(xlab)
plt.ylabel(ylab)
plt.title(xlab+' vs. '+ylab)
pp.savefig()
pp.close()
plt.show()
