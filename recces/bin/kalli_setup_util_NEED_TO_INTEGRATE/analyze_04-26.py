from data_new import *
import turner_util_new as util

#curr_wt = [0.73, 0.1, 0.0071, 0, 4.26, 2.46, 0.25, 0, 0.0, 4.54]
curr_wt = [0.73, 0.1, 0.0071, 0, 4.26, 2.46, 0.25, 0, 1.54, 4.54]

# Figure out the canonical NN free energies
seqs = [
    ('cc', 'gg')]

#seqs = [
#    ('aa', 'uu'), ('au', 'au'), ('ac', 'gu'), ('ag', 'cu'), ('ua', 'ua'),
#    ('uc', 'ga'), ('ug', 'ca'), ('cc', 'gg'), ('cg', 'cg'), ('gc', 'gc')]


# sim = SingleSimulation('ST/', curr_wt)


#def compute_freeE(DoS, scores, total_volume):
#    kT_in_kcal = 0.593
#    kT_list = np.arange(0.3, 2, 0.01)
#    n_kT = kT_list.shape[0]
#    boltzmann_factor = np.exp(-np.tile(scores, (n_kT, 1)).T / kT_list).T
#    DoS_volume = np.sum(DoS)
#
#    Z_list = total_volume * np.sum(DoS / DoS_volume * boltzmann_factor, axis=1)
#    A_list = -kT_in_kcal * np.log(Z_list)
#    E_list = (
#        np.sum(DoS * scores * boltzmann_factor, axis=1) /
#        np.sum(DoS * boltzmann_factor, axis=1) / kT_list * kT_in_kcal)
#    minus_TS_list = A_list - E_list
#    return kT_list, A_list, E_list, minus_TS_list

def seq_parse(seq):
    in_brac = False
    parsed_seq = []
    if not seq:
        return parsed_seq
    i = 0
    for j, c in enumerate(seq):
        if j == 0:
            continue
        if in_brac:
            if c == ']':
                in_brac = False
        elif c == '[':
            in_brac = True
        else:
            parsed_seq.append(seq[i:j])
            i = j
    parsed_seq.append(seq[i:])
    return parsed_seq

def get_delta_G(seq1, seq2 ):
#def get_delta_G(seq1, seq2, filename='freeE_orig.npy'):
    seq1_parsed = seq_parse(seq1)
    seq2_parsed = seq_parse(seq2)
    len1 = len(seq1_parsed)
    len2 = len(seq2_parsed)

    def freeE_from_seq(seq):
	# actually calculate this here
        sim = SingleSimulation('raw_data/%s/' %(seq), curr_wt)
        split_seq = seq.split('_')
        s1 = split_seq[0]
        s2 = ''
        if len(split_seq) > 1:
       	       s2 = split_seq[1]
        G = sim.value - np.log( util.torsion_volume(s1,s2) )
        return G
    #kT_list = np.load('raw_data/g%s/%s' % (seq1, filename))[:, 0]
    assert(0 < len1 <= 2 and 0 < len2 <= 2)
    if len1 != len2:  # Dangling end case
        dup1 = 'g%s_%sc' % (seq1, seq2)
        if len1 > len2:
            strand0 = ''.join(['g'] + seq1_parsed[:-1])
            strand1 = 'g' + seq1
            dup0 = '%s_%sc' % (strand0, seq2)
        else:
            strand0 = ''.join(seq2_parsed[1:] + ['c'])
            strand1 = seq2 + 'c'
            dup0 = 'g%s_%s' % (seq1, strand0)
        freeE = (
            freeE_from_seq(dup1) - freeE_from_seq(strand1) -
            freeE_from_seq(dup0) + freeE_from_seq(strand0))
    else:  # Duplex case
        def freeE_duplex():
            strand0 = ''.join(['g'] + seq1_parsed[:-1])
            strand1 = ''.join(seq2_parsed[1:] + ['c'])
            strand2 = 'g' + seq1
            strand3 = seq2 + 'c'
            dup0 = '%s_%s' % (strand0, strand1)
            dup1 = 'g%s_%sc' % (seq1, seq2)
            freeE = (
                freeE_from_seq(dup1) - freeE_from_seq(strand2) -
                freeE_from_seq(strand3) - freeE_from_seq(dup0) +
                freeE_from_seq(strand0) + freeE_from_seq(strand1))
            # for cooper 08/30
            print dup1
            print freeE_from_seq(dup1)
            print strand2
            print freeE_from_seq(strand2)
            print strand3
            print freeE_from_seq(strand3)
            print dup0
            print freeE_from_seq(dup0)
            print strand0
            print freeE_from_seq(strand0)
            print strand1
            print freeE_from_seq(strand1)
            return freeE
        freeE = freeE_duplex()
        print "freeE for ",seq1," ", seq2, ":", freeE

        seq1_parsed, seq2_parsed = seq2_parsed, seq1_parsed
        seq1 = ''.join(seq1_parsed)
        seq2 = ''.join(seq2_parsed)
        freeE_flipped = freeE_duplex()
        print "freeE for ",seq1," ", seq2, ":", freeE_flipped

        freeE += freeE_flipped
        freeE *= 0.5
    return freeE
    #return kT_list, freeE

for seq in seqs:
    print seq
    dG = get_delta_G( seq[0], seq[1])
    print KT_IN_KCAL*dG




# From Rhiju's script
##!/usr/bin/python2.7
#from recces.data import SingleSimulation, KT_IN_KCAL, N_SCORE_TERMS
#from recces.util import torsion_volume
#from rhiju_setup_util import get_seqs
#from scipy.misc import logsumexp
#from math import log
#from os import getcwd, chdir
#
## recces wts.
#curr_wt = [0.73, 0.1, 0.0071, 0, 4.26, 2.46, 0.25, 0, 1.54, 4.54]
#CWD = getcwd()
#val = {}
#for dir in ['cc/','gcc/','gg/','gg_cc/','ggc/','ggc_gcc/']:
#    tag = dir[:-1]
#    chdir( dir )
#    print getcwd()
#    sim = SingleSimulation( 'ST/', curr_wt )
#    seq1,seq2 = get_seqs( dir )
#    val[tag] = sim.value - log( torsion_volume( seq1, seq2 ) )
#    print tag, val[tag], val[tag] * KT_IN_KCAL
#    chdir( CWD )
#
#print 'ggc_gcc - ggc - gcc - (gg_cc - gg - cc )'
#val_final = val['ggc_gcc'] - val['ggc'] - val['gcc'] - ( val['gg_cc'] - val['gg'] - val['cc'] )
#print val_final, val_final * KT_IN_KCAL
#
##  LocalWords:  sim
