#!/usr/bin/python2.7
# From Rhiju's script
from recces.data import SingleSimulation, KT_IN_KCAL, N_SCORE_TERMS
from recces.util import torsion_volume
from rhiju_setup_util import get_seqs
from scipy.misc import logsumexp
from math import log
from os import getcwd, chdir

# recces wts.
curr_wt = [0.73, 0.1, 0.0071, 0, 4.26, 2.46, 0.25, 0, 1.54, 4.54]
CWD = getcwd()
val = {}
for dir in ['cc/','gcc/','gg/','gg_cc/','ggc/','ggc_gcc/']:
    tag = dir[:-1]
    chdir( dir )
    print getcwd()
    sim = SingleSimulation( 'ST/', curr_wt )
    seq1,seq2 = get_seqs( dir )
    val[tag] = sim.value - log( torsion_volume( seq1, seq2 ) )
    print tag, val[tag], val[tag] * KT_IN_KCAL
    chdir( CWD )

print 'ggc_gcc - ggc - gcc - (gg_cc - gg - cc )'
val_final = val['ggc_gcc'] - val['ggc'] - val['gcc'] - ( val['gg_cc'] - val['gg'] - val['cc'] )
print val_final, val_final * KT_IN_KCAL

#  LocalWords:  sim
