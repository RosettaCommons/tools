#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Loads experimental epitope information in a database for efficient lookup at design time.

@author: Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu 
"""

import argparse, collections, os, sys

from epilib.epitope_predictor import EpitopeScore, EpitopeMap
from epilib.epitope_predictor_matrix import EpitopePredictorMatrix, Propred
from epilib.epitope_csv import EpitopeCSV
from epilib.epitope_database import EpitopeDatabase
from epilib.expt_data import IEDBData
from epilib.netmhcII import NetMHCII

# TODO potential additional functionality:
# * merge datasets -- requires reconciling conflicting data, and note that currently only positive hits are stored

def get_core_hits(scores, pred, option):
    """Uses epitope predictor to generate 9mer core hits for the peptide epitope scores.
    scores: peptide => allele => score
    pred: EpitopePredictor returning 9mer
    option: 'all', predicted_good', 'predicted_best' -- which core(s) #TODO enum?
    returns: core => [allele] such that there is evidence that the core is for a peptide binder against allele
    """
    # TODO: currently just indicates whether there's evidence for a core binding an allele; nothing relative about strength of evidence
    # if want to do something quantitative, base it on scores[peptide][allele]
    core_hits = collections.defaultdict(set) # core => [allele]
    
    for peptide in scores:
        if pred:
            pep_epimap = pred.score_protein(peptide)
        for allele in scores[peptide]:
            if option == 'all':
                for i in range(len(peptide)-9): 
                    core_hits[peptide[i:i+9]].add(allele)
            elif option == 'predicted_good':
                a = pred.alleles.index(allele) # TODO: ugh -- details is just a list in order of alleles
                ncores = 0
                for core in pep_epimap.peptides:
                    if pep_epimap.peptide_score(core).details[a] is not None: # TODO: Propred does this, but not NetMHC (it also doesn't store 9mers, so fix simultaneously)
                        core_hits[core].add(allele)
                        ncores += 1
                # print(peptide,allele,ncores,'cores')
            else:
                # TODO: implement
                raise Exception('not implemented')
    #print(len(core_hits),'cores')
    return core_hits

def setup_parser():
    """Creates an ArgumentParser and sets up all its arguments.
    returns: ArgumentParser"""
    parser = argparse.ArgumentParser(description='Establishes a database of precomputed epitope scores, to enable efficient lookup at design time', 
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog="""Additional notes:
- Allele subsetting (--allele_set / --alleles) is required for mysql and allowed for csv import
- Allele naming is wrt whatever is used in the input file/db (currently IEDB)
- IEDB includes names that have only DRB as well as those paired up with DRA (e.g., HLA-DRB1*01:01 and HLA-DRA*01:01/DRB1*01:01); the predefined sets include both
- IEDB also includes some not completely specified names (e.g., HLA-DR1); these are ignored but presumably could be pulled in with specific naming
- If the specified "--db" or "--csv" exists, an error is raised unless "--overwrite" or "--augment" is indicated (augment currently only supported for db)
- If an epitope predictor is used, an attempt will be made to convert allele names to those used by the predictor
- NetMHCII identifies one core per peptide (may be of low confidence), but matrices can identify multiple (or none)
    """)

    # where to get the data
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--iedb_csv', help='name of IEDB-format csv export')
    source.add_argument('--iedb_mysql', help='name of mysql database holding IEDB-format tables')
    source.add_argument('--iedb_fresh_mysql', help='name of mysql database into which IEDB-format tables will be downloaded (via curl) and stored, to enable processing (note: blows away existing IEDB tables!)')

    parser.add_argument('--mysql_user', help='username for connecting to mysql database (default: %(default)s)', default='root')
    parser.add_argument('--mysql_pw', help='password for connecting to mysql database (default: %(default)s)', default=None)

    # from peptide to consistuent core 9mer(s)
    parser.add_argument('--cores', help='which 9mer core(s) of a peptide to include (default: %(default)s)', choices=['all','predicted_good','predicted_best'], default='all')
    # subset of alleles to take from data
    alleles = parser.add_mutually_exclusive_group()
    alleles.add_argument('--allele_set', help='name of predefined set of alleles', choices=['test', 'greenbaum11', 'paul15'])
    alleles.add_argument('--alleles', help='comma-separated list of allele names')

    # epitope predictor
    pred = parser.add_mutually_exclusive_group()
    pred.add_argument('--matrix', help='generic epitope predictor matrix filename')
    #pred.add_argument('--netmhcii', action='store_true', help='use netmhcII executable')
    pred.add_argument('--pred_csv', help='name of csv from which to load epitope predictions') # TODO: different name from score.py since output is also db -- confusing?
    pred.add_argument('--pred_db', help='name of database from which to load epitope predictions') # TODO: different name from score.py since output is also db -- confusing?
    pred.add_argument('--propred', action='store_true', help='use propred matrices')
    # epitope predictor parameters
    parser.add_argument('--epi_thresh', help='epitope predictor threshold (default: %(default).2f)', type=int, default=5)
    parser.add_argument('--noncanon', help='how to treat letters other than the 20 canonical AAs (default: %(default)s)', choices=['error', 'silent', 'warn'])
    parser.add_argument('--netmhcii_score', help='type of score to compute (default %(default)s)', choices=['rank','absolute'], default='rank')

    # the database
    out = parser.add_mutually_exclusive_group()
    out.add_argument('--db', help='name of sqlite3 database to store epitope information (create or augment)')
    out.add_argument('--csv', help='name of csv file to store epitope information (create)')
    oa = parser.add_mutually_exclusive_group()
    oa.add_argument('--augment', help='augment existing db (not currently supported for csv)', action='store_true')
    oa.add_argument('--overwrite', help='overwrite existing db or csv', action='store_true')

    return parser

def main(args, argv):
    """Sets up a data-based epitope database based on the parsed args. (Say that 5 times fast.)
    args: argparse.Namespace (to use)
    argv: original sys.argv (to store)"""
    
    # load data
    data = None
    data_allele_list = None if args.alleles is None else args.alleles.split(',')
    if args.iedb_csv is not None:
        data = IEDBData.from_csv(args.iedb_csv, data_allele_list, args.allele_set)
    elif args.iedb_mysql is not None:
        data = IEDBData.from_mysql(args.iedb_mysql, data_allele_list, args.allele_set, args.mysql_user, args.mysql_pw)
    elif args.iedb_fresh_mysql is not None:
        IEDBData.download_mysql(args.iedb_fresh_mysql, user=args.mysql_user, pw=args.mysql_pw)
        if data_allele_list is None and args.allele_set is None:
            # only goal was to download the database, so bail
            return
        # else continue with the fresh db
        data = IEDBData.from_mysql(args.iedb_fresh_mysql, data_allele_list, args.allele_set, user=args.mysql_user, pw=args.mysql_pw)
    
    # alleles
    alleles = list(data.std_alleles)
    
    # epitope predictor
    pred = None # default is not to use prediction
    if args.propred:
        pred = Propred.load()
    elif args.matrix is not None:
        pred = EpitopePredictorMatrix.load(args.matrix)
#    elif args.netmhcii:
#        pred = NetMHCII(score_type=args.netmhcii_score[0])
#        raise Exception('not implemented yet -- need to extract cores')
    elif args.pred_csv is not None:
        pred = EpitopeCSV.for_reading(args.pred_csv)
    elif args.pred_db is not None:
        pred = EpitopeDatabase.for_reading(args.pred_db)
    if pred is not None:
        pred.thresh = args.epi_thresh
        if pred.peptide_length != 9:
            raise Exception('need to predict 9mer cores')
        if args.allele_set is not None:
            if args.allele_set not in pred.allele_sets: 
                raise Exception('allele_set '+args.allele_set+' not supported')
            pred.filter_alleles(pred.allele_sets[args.allele_set])
        elif args.alleles is not None:
            # TODO: convert allele names from those for IEDB to those for predictor, by way of data.std_alleles and corresponding std_name in predictor
            # a bit gross -- any better way?
            raise Exception('not implemented')

    # TODO: when extending to NetMHC, at this point should batch predict for the identified peptides
    # will then need to extract predicted core
    # could also do a two-stage thing as with gen_db -- save out a file and let user run separately, load into db, and continue with that db
    
    # which cores to store
    if args.cores == 'predicted_good' or args.cores == 'predicted_best':
        if pred is None:
            raise Exception('if going to use predicted cores (--cores predicted_*), must specify an epitope predictor')
        # TODO: NetMHC only returns one predicted core, so disallow good or at least print warning
    core_hits = get_core_hits(data.scores, pred, args.cores)
    
    # make scores for the cores
    cores = list(core_hits.keys())
    # TODO: currently a binary score (is there a hit or not); could scale based on amount of evidence, # alleles, whatever
    scores = [EpitopeScore(1, [1 if a in hits else None for a in alleles]) for hits in core_hits.values()]
    epimap = EpitopeMap(9, alleles, cores, scores)
        
    # TODO: seems like it would be nice to have a pointer from the core back to it source, from which we could then get details of the data
    
    out = None
    if args.db is not None:
        if os.path.exists(args.db):
            if args.overwrite:
                os.remove(args.db)
            elif not args.augment:
                raise Exception('database '+args.db+' already exists; specify --augment or --overwrite')
        out = EpitopeDatabase.for_writing(args.db, ' '.join(argv), alleles, 9)
    elif args.csv is not None:
        if os.path.exists(args.csv):
            if args.overwrite:
                os.remove(args.csv)
            else:
                raise Exception('csv file '+args.csv+' already exists; specify --overwrite')
        out = EpitopeCSV.for_writing(args.csv, alleles, 9)
    if out is None:
        print('no output generated')
    else:
        out.save_scores(epimap)
        
if __name__ == '__main__':
    parser = setup_parser()
    args = parser.parse_args()
    # main(args, sys.argv)
    try:
        main(args, sys.argv)
    except Exception as e:
        print('ERROR', e)
        print()
        parser.print_usage()
