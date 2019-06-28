#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute background distributions for YAEPII models, in order to convert predicted values to percentile ranks

uniprot-reviewed%3Ayes+AND+proteome%3Aup000005640.fasta

@author: Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu 
"""

import argparse, os
from epilib.yaepII import YAEPII
from epilib.epitope_predictor import EpitopeScore, EpitopeMap
from epilib.epitope_database import EpitopeDatabase
from epilib.sequence import load_pep, load_fsa

def setup_parser():
    """Creates an ArgumentParser and sets up all its arguments.
    returns: ArgumentParser"""

    parser = argparse.ArgumentParser(description='Computes background distributions for YAEPII models', 
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog="""Additional notes:
    - if use a db, everything in the db contributes to the dist
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
    parser.add_argument('--dist_file', help='filename (no directory) for resulting distribution (default: %(default))', default='dist.txt')
    parser.add_argument('--overwrite', help='must set this in order to overwrite existing distribution file', action='store_true')

    return parser

def compute_one_allele(args, seqs, allele):
    """Handles calculation for one specified allele.
    args: ArgumentParser""" 

    # directory for this allele   
    base = args.model_base if args.model_base is not None else 'yaepII'
    allele_dir = base+'/'+allele+'/'

    # check that out file is legit    
    fn = allele_dir + args.dist_file
    if os.path.exists(fn):
        if args.overwrite:
            print('overwriting',fn)
        else:
            raise Exception('distribution "'+fn+'" already exists; specify --overwrite if you want to overwrite it')

    model = YAEPII.load_frozen([allele], args.model_base).alleles[0] # TODO: little awkward since doing an allele at a time, but ok?
    
    # database
    if args.build_db:
        db = EpitopeDatabase.for_writing(allele_dir + '/' + args.db, 'yaepII_bg', alleles=[allele], peptide_length=15)
    else:
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
    # TODO: don't yet know what format will be best to use in Rosetta, so just raw sorted scores
    dist = sorted(scores)

    # save
    with open(fn, 'w') as out:
        print(dist, file=out)

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
