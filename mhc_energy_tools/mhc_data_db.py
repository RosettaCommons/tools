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
#from epilib.netmhcII import NetMHCII
from epilib.ext_citations import CitationTracker

# TODO potential additional functionality:
# * general control over allele synonyms (w/o allele file)
# * the epitope prediction stuff is quite skeletal, and doesn't yet support NetMHCII
# * finer control over MHC assays
# * beyond MHC data to T cell data (in a different table in IEDB)
# * seems like it could be nice to have a pointer from the stored peptide back to the source data, for subsequent analysis   
# * scores are currently binary (is there a hit or not); could scale based on amount of evidence
# * merge datasets -- requires reconciling conflicting data, and note that currently only positive hits are stored

def get_core_hits(expt, pred, option):
    """Uses epitope predictor to generate 9mer core hits for the peptide epitope scores.
    expt: peptide => allele => score
    pred: EpitopePredictor returning 9mer
    option: 'all', predicted_good', 'predicted_best', 'full' -- which core(s) #TODO enum?
    returns: core => [allele] such that there is evidence that the core is for a peptide binder against allele
    """
    # TODO: currently just indicates whether there's evidence for a core binding an allele; nothing relative about strength of evidence
    # if want to do something quantitative, base it on scores[peptide][allele]
    core_hits = collections.defaultdict(set) # core => [allele]
    
    for peptide in expt:
        if pred:
            pep_epimap = pred.score_protein(peptide)
        for allele in expt[peptide]:
            if option == 'all':
                for i in range(len(peptide)-9): 
                    core_hits[peptide[i:i+9]].add(allele)
            elif option == 'predicted_good':
                # prediction for experimentally identified allele
                a = pred.alleles.index(allele) # TODO: ugh -- details is just a list in order of alleles
                ncores = 0
                for core in pep_epimap.peptides:
                    if pep_epimap.peptide_score(core).details[a] is not None: # TODO: Propred does this, but not NetMHC (it also doesn't store 9mers, so fix simultaneously)
                        core_hits[core].add(allele)
                        ncores += 1
                # print(peptide,allele,ncores,'cores')
            elif option == 'full':
                core_hits[peptide].add(allele)
            else:
                # TODO: implement
                raise Exception(option+' not implemented')
    #print(len(core_hits),'cores')
    return core_hits

def setup_parser():
    """Creates an ArgumentParser and sets up all its arguments.
    returns: ArgumentParser"""
    parser = argparse.ArgumentParser(description='Establishes a database of precomputed epitope scores, to enable efficient lookup at design time', 
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog="""Additional notes:
- Allele subsetting (--allele_set / --alleles / --allele_file) is required for mysql and allowed for csv import.
- The allele_set 'hlaII' indicates all alleles known to be of HLA class II. For mysql, that's done by query; for csv, using the provided file 'data/hlaII.csv' which was exported from IEDB but may not be the latest.
- Allele naming is wrt whatever is used in the input file/db (currently IEDB)
- IEDB includes names that have only DRB as well as those paired up with DRA (e.g., HLA-DRB1*01:01 and HLA-DRA*01:01/DRB1*01:01); the predefined sets include both.
- IEDB also includes some not completely specified names (e.g., HLA-DR1); these are ignored but presumably could be pulled in with specific naming
- The allele_file allows precise control over which alleles to use. Alternatively, the IEDB web server can be used to pull down data for just the alleles of interest.
- Allele_file (e.g., 'data/iedb-hlaII.csv'): one allele name per row. Header row required; fixed column format anyway so names don't matter. Optional second column (comma separated) allows specifying of synonyms to be merged to a common name; e.g., HLA-DRA*01:01/DRB*01:01 to simply HLA-DRB*01:01, which is a different row). Optional third column gives the name as used by an epitope predictor
- If the specified "--db" or "--csv" exists, an error is raised unless "--overwrite" or "--augment" is indicated (augment currently only supported for db)
- Peptides are of variable length, so are converted into one or more 9mer core(s) considered risky -- all of them, or just the best or "good enough" ones according to an epitope predictor. The 'full' option saves out the whole peptide without 'core'ing it; this can't be used in the current mhc_epitope score but might be useful.
- If an epitope predictor is used, only the experimentally-evaluated allele will be predicted, and an attempt will be made to convert allele names to those used by the predictor. The allele_file can specify such a mapping.
- NetMHCII (currently not supported) identifies one core per peptide (may be of low confidence), but matrices can identify multiple (or none)
- Assay restrictions are currently only supported for queries to the local database, and currently are either just all/none for binding and elution separately. This follows IEDB's breakdown. The IEDB website allows finer-grained filtering; this script could be extended to do likewise.
    """)

    # where to get the data
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--iedb_csv', help='name of IEDB-format csv export')
    source.add_argument('--iedb_mysql', help='name of mysql database holding IEDB-format tables')
    source.add_argument('--iedb_fresh_mysql', help='name of mysql database into which IEDB-format tables will be downloaded (via curl) and stored, to enable processing (note: blows away existing IEDB tables!)')
    parser.add_argument('--mysql_user', help='username for connecting to mysql database (default: %(default)s)', default='root')
    parser.add_argument('--mysql_pw', help='password for connecting to mysql database (default: %(default)s)', default=None)

    # types of assay (following IEDB's breakdown)
    parser.add_argument('--assay_mhc_ligand_binding', help='which forms of MHC binding data to use; only supported with mysql interface (default: %(default)s)', choices=['all','none'], default='all') # TODO: break down subtypes
    parser.add_argument('--assay_mhc_ligand_elution', help='which forms of MHC ligand elution binding data to use; only supported with mysql interface (default: %(default)s)', choices=['all','none'], default='all') # TODO: break down subtypes

    # from peptide to consistuent core 9mer(s)
    parser.add_argument('--cores', help='which 9mer core(s) of a peptide to include (default: %(default)s).  all means all 9mers, predicted_good/_best selects good binders/best binder core(s), and full leaves the intact peptide (does not reduce to a 9mer).', choices=['all','predicted_good','predicted_best','full'], default='all')

    # subset of alleles to take from data
    alleles = parser.add_mutually_exclusive_group()
    alleles.add_argument('--allele_set', help='name of predefined set of alleles', choices=['test', 'greenbaum11', 'paul15', 'southwood98', 'hlaII'])
    alleles.add_argument('--alleles', help='comma-separated list of allele names')
    alleles.add_argument('--allele_file', help='name of csv-file with allele names, format in "additional notes" section')

    # epitope predictor
    pred = parser.add_mutually_exclusive_group()
    pred.add_argument('--matrix', help='generic epitope predictor matrix filename')
    #pred.add_argument('--netmhcii', action='store_true', help='use netmhcII executable')
    pred.add_argument('--pred_csv', help='name of csv from which to load epitope predictions')
    pred.add_argument('--pred_db', help='name of database from which to load epitope predictions')
    pred.add_argument('--propred', action='store_true', help='use propred matrices')
    # epitope predictor parameters
    parser.add_argument('--epi_thresh', help='epitope predictor threshold (default: %(default).2f)', type=int, default=5)
    parser.add_argument('--noncanon', help='how to treat letters other than the 20 canonical AAs (default: %(default)s)', choices=['error', 'silent', 'warn'])
    #parser.add_argument('--netmhcii_score', help='type of score to compute (default %(default)s)', choices=['rank','absolute'], default='rank')

    # output db / file
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
    
    # Citation manager object to track what features are being used
    cite_manager = CitationTracker()
    # load data
    data = None
    data_allele_list = None if args.alleles is None else args.alleles.split(',')
    if args.iedb_csv is not None:
        cite_manager.add_citation(cite_manager.iedb)
        if args.assay_mhc_ligand_binding != 'all' or args.assay_mhc_ligand_elution != 'all': print('csv currently cannot restrict assay type; ignoring that flag')
        data = IEDBData.from_csv(args.iedb_csv, data_allele_list, args.allele_set)
    elif args.iedb_mysql is not None:
        cite_manager.add_citation(cite_manager.iedb)
        data = IEDBData.from_mysql(args.iedb_mysql, args.assay_mhc_ligand_binding, args.assay_mhc_ligand_elution, data_allele_list, args.allele_set, user=args.mysql_user, pw=args.mysql_pw)
    elif args.iedb_fresh_mysql is not None:
        cite_manager.add_citation(cite_manager.iedb)
        IEDBData.download_mysql(args.iedb_fresh_mysql, user=args.mysql_user, pw=args.mysql_pw)
        if data_allele_list is None and args.allele_set is None:
            # only goal was to download the database, so bail
            cite_manager.output_citations()
            return
        # else continue with the fresh db
        data = IEDBData.from_mysql(args.iedb_fresh_mysql, args.assay_mhc_ligand_binding, args.assay_mhc_ligand_elution, data_allele_list, args.allele_set, user=args.mysql_user, pw=args.mysql_pw)
    
    # alleles
    std_alleles = sorted(set(data.std_alleles.values()))
    
    # epitope predictor
    pred = None # default is not to use prediction
    if args.propred:
        pred = Propred.load()
        cite_manager.add_citation(cite_manager.propred)
    elif args.matrix is not None:
        #We don't have any generic matrix strategy.  If added, add the citation.
        pred = EpitopePredictorMatrix.load(args.matrix)
#    elif args.netmhcii:
#        cite_manager.add_citation(cite_manager.netmhcii)
#        pred = NetMHCII(score_type=args.netmhcii_score[0])
#        raise Exception('not implemented yet -- need to extract cores')
    elif args.pred_csv is not None:
        #Generic, so we can't add a citation
        pred = EpitopeCSV.for_reading(args.pred_csv)
    elif args.pred_db is not None:
        #Generic, so we can't add a citation
        pred = EpitopeDatabase.for_reading(args.pred_db)
    if pred is not None:
        pred.thresh = args.epi_thresh
        if pred.peptide_length != 9:
            raise Exception('need to predict 9mer cores')
        if args.allele_set is not None:
            if args.allele_set not in pred.allele_sets: 
                # TODO: if predictor is a CSV/DB, it doesn't have an allele_set per se.  This becomes a bit esoteric, though, since you are unlikely to make a db to then pick cores.
                raise Exception('allele_set '+args.allele_set+' is not supported by ' + pred.name + ', the Predictor you are using for picking cores.  You must select one of ' + pred.name + "'s supported allele_sets: " + str(list(pred.allele_sets.keys())))
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
    scores = [EpitopeScore(1, [1 if a in hits else None for a in std_alleles]) for hits in core_hits.values()]
    epimap = EpitopeMap(9, std_alleles, cores, scores)
        
    out = None
    if args.db is not None:
        if os.path.exists(args.db):
            if args.overwrite:
                os.remove(args.db)
            elif not args.augment:
                raise Exception('database '+args.db+' already exists; specify --augment or --overwrite')
        out = EpitopeDatabase.for_writing(args.db, ' '.join(argv), std_alleles, 9)
    elif args.csv is not None:
        if os.path.exists(args.csv):
            if args.overwrite:
                os.remove(args.csv)
            else:
                raise Exception('csv file '+args.csv+' already exists; specify --overwrite')
        out = EpitopeCSV.for_writing(args.csv, std_alleles, 9)
    if out is None:
        print('no output generated')
    else:
        out.save_scores(epimap)
    
    cite_manager.output_citations()

if __name__ == '__main__':
    parser = setup_parser()
    args = parser.parse_args()
    try:
        main(args, sys.argv)
    except Exception as e:
        print('ERROR', e)
        print()
        parser.print_usage()
        # raise(e)
