#!/usr/bin/env python

from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse

parser = argparse.ArgumentParser(description='Plot scores terms from silent files.')
parser.add_argument('silent_file', help='Name of the Silent File')
parser.add_argument('-tag', help='Tag of Output File', default='score_vs_rmsd')
parser.add_argument('-x', help='X-axis variable', default='NAT_rmsd')
parser.add_argument('-y', help='Y-axis variable', default='score')
parser.add_argument('-nstruct', default=None)
parser.add_argument('-cycles', default=None)
parser.add_argument('-fit',default=True)
args=parser.parse_args()

silent_file=args.silent_file
tag=args.tag
fit=args.fit

if args.cycles:
	tag=args.cycles+'_cycles_'+tag
if args.nstruct:
	tag=args.nstruct+'_models_'+tag
if '.out' in silent_file:
	silent_tag=silent_file.replace('.out','_')
if '_rebuild' in silent_tag:
	silent_tag=silent_tag.replace('_rebuild','')
tag=silent_tag+tag


out_pdf=tag+'.pdf'
out_png=tag+'.png'
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
details=None
models=None
cycles=None
if args.nstruct:
	models = args.nstruct+' models'
if args.cycles:
	cycles = args.cycles+' cycles'

if models and cycles:
	details='\n('+models+', '+cycles+')'
elif models:
	details='\n('+models+')'
elif cycles:
	details='\n('+cycles+')'
else:
	details=''


pp=PdfPages(out_pdf)
plot = plt.plot(xdata,ydata, marker='o', color='r', linestyle=' ')
plt.title(ylab+' vs. '+xlab+details)
plt.xlabel(xlab)
plt.ylabel(ylab)
if not fit:
	plt.xlim(0.0,9.0)
	plt.ylim(-40,0)
pp.savefig()
pp.close()
plt.show()
plt.savefig(out_png)
