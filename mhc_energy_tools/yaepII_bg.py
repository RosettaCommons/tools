#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute percentile ranks from background distributions for YAEPII models
The output files specify, for sampled score values (implicit, at steps of 0.001) the corresponding percentile ranks
Percentile ranks are in the "lower is better" order, i.e., 5 means a higher probability that 95% of the background peptides

uniprot-reviewed%3Ayes+AND+proteome%3Aup000005640.fasta

@author: Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu 
"""

import argparse, os
from epilib.yaepII import YAEPIIAlleleModel, YAEPII
from epilib.epitope_predictor import EpitopeScore, EpitopeMap
from epilib.epitope_database import EpitopeDatabase
from epilib.sequence import load_pep, load_fsa

import numpy as np
import scipy
import statsmodels.distributions.empirical_distribution as edf

def setup_parser():
    """Creates an ArgumentParser and sets up all its arguments.
    returns: ArgumentParser"""

    parser = argparse.ArgumentParser(description='Computes percentile ranks from background distributions for YAEPII models', 
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog="""Additional notes:
    - all files are stored within a model's directory
    - the distribution goes into ranks.txt, with one line for each score, stepping by 0.0001, giving the percentile rank (lower is better) for the score (e.g., line 1 is the percentile rank of score 0.001)
    - if using a db, augment it, and then everything in the db contributes to the calculation
                                         """)
    # where to get the sequences
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--fsa', help="name of file with one or more sequences in fasta-ish format, each including '>' line")
    source.add_argument('--pep', help="name of file with one sequence per line")

    # which allele(s)
    alleles = parser.add_mutually_exclusive_group(required=True)
    alleles.add_argument('--allele_set', help='name of predefined set of alleles', choices=YAEPII.allele_sets.keys())
    alleles.add_argument('--alleles', help='comma-separated list of allele names')

    # where to find models and store background distributions    
    parser.add_argument('--model_base', help='base directory for models. etc. (nested within) (default: %(default))', default='models')

    # outputs
    parser.add_argument('--build_db', help='create/augment a database with peptide scores, so build up over runs', action='store_true')
    parser.add_argument('--db', help='filename (no directory) of db to be created/augmented if build_db is true (default: %(default))', default='scores.db')
    parser.add_argument('--plot_dist', help='if set, plot score distribution for inspection (matplotlib required), in dist.pdf', action='store_true')
    parser.add_argument('--overwrite', help='must set this in order to overwrite existing file', action='store_true')

    return parser

def compute_one_allele(args, seqs, allele):
    """Handles calculation for one specified allele.
    args: ArgumentParser""" 

    # directory for this allele   
    base = args.model_base if args.model_base is not None else 'yaepII'
    allele_dir = base+'/'+allele+'/'

    # check that out file is legit    
    fn = allele_dir + 'ranks.txt'
    if os.path.exists(fn):
        if args.overwrite:
            print('overwriting',fn)
        else:
            raise Exception('ranks file "'+fn+'" already exists; specify --overwrite if you want to overwrite it')

    model = YAEPIIAlleleModel.load_frozen(allele, args.model_base, False)
    
    if args.build_db:
        # will augment database, then get all scores
        db = EpitopeDatabase.for_writing(allele_dir + '/' + args.db, 'yaepII_bg', alleles=[allele], peptide_length=15)
    else:
        # just current scores
        scores = []

    # scan peptides in sequences
    for seq in seqs:
        # include N- and C-terminal peptides padded to length 15
        padded = '---' + seq.seq + '---'
        peps = [padded[i:i+15] for i in range(len(padded)-14)]
        # TODO: batch score; if db, only for new ones
        pep_scores = [model.score_peptide(pep) for pep in peps]
        if args.build_db:
            db.save_scores(EpitopeMap(15, [allele], peps, [EpitopeScore(float(score), [float(score)]) for score in pep_scores]))
        else:
            scores.extend(pep_scores)       

    # get everything from db (including previous runs, if any)
    if args.build_db:
        scores = [s.value for s in db.load_scores().scores]
            
    # compute distribution
    # inspired by https://stackoverflow.com/questions/44132543/python-inverse-empirical-cumulative-distribution-function-ecdf
    scores_edf = edf.ECDF(scores)
    slope_changes = sorted(set(scores) | set([0,1]))
    scores_edf_values_at_slope_changes = [scores_edf(score) for score in slope_changes]
    interp_edf = scipy.interpolate.interp1d(slope_changes, scores_edf_values_at_slope_changes)
    sample_scores = np.arange(0, 1.001, 0.001)
    sample_pcts = (1-interp_edf(sample_scores))*100

    # save
    with open(fn, 'w') as out:
        for pct in sample_pcts: print(pct, file=out)

    if args.plot_dist:
        try:
            import matplotlib.pyplot as plt
        except:
            print('matplotlib is required for plotting, and could not be imported.  Make sure it is correctly installed.\n')
            exit()
        plt.clf()
        plt.hist(scores, color='lightgray')
        ax2 = plt.gca().twinx()
        ax2.plot(sample_scores, sample_pcts, 'b.')
        plt.savefig(fn[:-3]+'pdf')

def main(args):
    # sequences
    if args.fsa is not None:
        seqs = load_fsa(args.fsa)
    elif args.pep is not None:
        seqs = load_pep(args.pep)

    # alleles
    if args.allele_set is not None:
        # note that the allowed choices were limited to the dictionary keys, so should be legit
        alleles = YAEPII.allele_sets[args.allele_set]
    elif args.alleles is not None:
        alleles = args.alleles.split(',')
    if len(alleles) == 0:
        print('!!! no alleles specified -- this will be easy :)')
    
    # TODO: distribute across processors?
    for allele in alleles:
        print('=================',allele)
        compute_one_allele(args, seqs, allele)
    
# python yaepII_bg.py --fsa ~/Downloads/test.fsa --alleles DRB1_0101
# python yaepII_bg.py --fsa ~/Downloads/uniprof_human_reviewed.fsa --alleles DRB1_0101
if __name__ == '__main__':
    parser = setup_parser()
    args = parser.parse_args()
    main(args)
