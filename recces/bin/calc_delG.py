#!/usr/bin/env python

from recces.data import SingleSimulation,  N_SCORE_TERMS
from recces.util import KT_IN_KCAL
from math import log,exp
from os.path import dirname, abspath
from sys import argv
import re
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='run delG calculation via RECCES')
parser.add_argument('-seq', help='Sequence of the helix, in the format: gg_cc' )
parser.add_argument('-incorrect_wt_order', help='Tried new order of turner weights for rb_recces.',action='store_true', default=False )
parser.add_argument('-no_kcal_mol', help='Put out numbers in kT, not kcal/mol',action='store_true', default=False )
parser.add_argument('-weight_sets',help='File with other weight sets to test' )
args = parser.parse_args()

# turner.wts was apparently in this order (see kalli dir 04-26)
curr_wt = [0.73, 0.1, 0.0071, 0, 4.26, 2.46, 0.25, 0, 1.54, 4.54]
# in early rb_recces tests, I ended up trying this order, but I think it was wrong. -- rhiju, december 2016
if args.incorrect_wt_order: curr_wt = [0.73, 0.1, 0.0071, 0.25, 4.54, 4.26, 0.0001, 1.54, 0.0001, 2.46 ]

sim = SingleSimulation('./', curr_wt, name=args.seq)
val = sim.value
if not args.no_kcal_mol: val *= KT_IN_KCAL

print "RECCES deltaG ",
if not args.no_kcal_mol: print "(kcal/mol)",
print ':',val
#if ( sim.value < 0.0 ):
#    print "with -1 correction, kcal/mol:", -KT_IN_KCAL * log( exp( -sim.value ) - 1 )
print


if ( args.weight_sets != None ):
    # WARNING! Not tested.
    # following function from cooper:
    # Establish the list of weights to be applied:
    f = open( dirname( abspath( argv[0] ) ) + '/weight_sets.txt','r+')
    A = []
    for line in f:
        A+=[line]
    weight_list=[]
    for strang in A:
        weight_list+=[re.findall(r"[-+]?\d*\.\d+|\d+",strang)]
    f.close

    # WARNING! Not tested.
    # from cooper
    # curr_wt = [0.73, 0.1, 0.0071, 0, 4.26, 2.46, 0.25, 0, 1.54, 4.54]
    reorder = [ 1, 2, 3, 7, 10, 5, 4, 9, 8, 6] # verified by looking at score_type output list and score_type order in turner.wts
    vals = []
    print "Computing average over other weights..."
    for weight in weight_list[0:1000:20]:
        if len(weight) == 11:
            sim.reweight( [float(weight[i]) for i in reorder ] )
            print "deltaG(corr): ", sim.value * KT_IN_KCAL
            vals.append( sim.value * KT_IN_KCAL)
        else:
            continue

    print
    print "deltaG(corr): ", np.mean( np.array( vals ) ), "+/-", np.std( np.array( vals ) ), " from ", len( vals ), " weight sets"
