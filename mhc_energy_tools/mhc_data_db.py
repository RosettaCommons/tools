#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Loads experimental epitope information in a database for efficient lookup at design time.

@author: Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu 
"""

import argparse, collections, csv

from epilib.epitope_predictor import EpitopeScore, EpitopeMap
from epilib.epitope_predictor_matrix import EpitopePredictorMatrix, Propred
from epilib.epitope_database import EpitopeDatabase
from epilib.netmhcII import NetMHCII
from epilib.sequence import AAs

# TODO added functionality
# merge datasets -- requires reconciling conflicting data, and note that currently only positive hits are stored

def internal_allele_name(iedb_name):
    # TODO: gross hack for now
    # need to further standardize allele names globally; see also started NetMHCII.std_name
   return iedb_name \
    .replace('HLA-','') \
    .translate(str.maketrans({'/':'_', '-':'_', '*':'_', ':':None})) \
    .replace('DRA_0101_','') # TODO: just noticed then one; other DRA alleles?

def load_iedb_csv(filename, iedb_alleles):
    measurements = {} # peptide => allele => [qual_measure]  TODO: keep more info?
    with open(filename) as infile:
        # IEDB splits the header over two lines, so concatenate them as use them as the field names
        head1 = infile.readline().strip().split(',')
        head2 = infile.readline().strip().split(',')
        nrows = 0
        for row in csv.DictReader(infile, fieldnames=[h1+' '+h2 for h1,h2 in zip(head1,head2)]):
            nrows += 1
            peptide = row['Epitope Description']
            iedb_allele = row['MHC Allele Name']
            if iedb_allele not in iedb_alleles:
                print('skipping allele',iedb_allele)
                continue
            internal_allele = internal_allele_name(iedb_allele)
            # only deal with peptides containing at least 9 amino acids
            if len(peptide)<9 or any(aa not in AAs for aa in peptide): continue
            qual_measure = row['Assay Qualitative Measure']
            if peptide not in measurements:
                measurements[peptide] = {internal_allele:[qual_measure]}
            elif internal_allele not in measurements[peptide]:
                measurements[peptide][internal_allele] = [qual_measure]
            else:
                measurements[peptide][internal_allele].append(qual_measure)
    print(nrows,'rows','=>',len(measurements),'peptides')
    return measurements

"""Notes about IEDB mysql; move into wiki docs when settled

* on my mac, installed mariadb via homebrew

> mysql.server start

* downloaded iedb sql  http://www.iedb.org/database_export_v3.php => http://www.iedb.org/downloader.php?file_name=doc/iedb_public.sql.gz

* imported that into a database
> mysql -uroot -e "create database iedb;"
> mysql -uroot iedb < iedb_public.sql

* installed one of the python packages to connect to the db
> pip install mysql-connector-python

* example query https://help.iedb.org/hc/en-us/articles/114094146451-Select-all-MHC-binding-assays-for-class-II-human-epitopes-with-a-SQL-queryAP

select m.mhc_bind_id into outfile "~/Downloads/bind.csv" from mhc_bind m where m.mhc_allele_name in ('HLA-DRB1*01:01', 'HLA-DRA*01:01/DRB1*01:01') order by m.mhc_bind_id;
select m.mhc_elution_id into outfile "~/Downloads/elute.csv" from mhc_elution m where m.mhc_allele_name in ('HLA-DRB1*01:01','HLA-DRA*01:01/DRB1*01:01') order by m.mhc_elution_id;
-> combining these gives almost the same as what I get via web query (just restricting HLA type to 0101); web query has two additional sets of epitopes, apparently from two additional publications (recently added?)

* fields http://curationwiki.iedb.org/wiki/index.php/Data_Field_Descriptions
"""

def load_iedb_mysql(dbname, iedb_alleles):
    # TODO: currently just binding data; also get elution from separate table
    # TODO: t cell data?
    # TODO: allow filtering on assay details
    query = ('select e.linear_peptide_seq as peptide, m.mhc_allele_name as allele, m.as_char_value as qual_measure '
             'from mhc_bind m, curated_epitope ce, epitope_object eo, epitope e '
             'where m.curated_epitope_id=ce.curated_epitope_id and ce.e_object_id = eo.object_id and eo.epitope_id=e.epitope_id '
             'and length(e.linear_peptide_seq)>=9 '
             'and m.mhc_allele_name in ') + '(' + ','.join('"'+a+'"' for a in iedb_alleles) + ')'
    print(query)
    
    import mysql.connector
    connection = mysql.connector.connect(database=dbname, user='root') # TODO: custom username and password
    cursor = connection.cursor()
    cursor.execute(query)
    # TODO: error handling
    
    measurements = {} # peptide => allele => [qual_measure]  TODO: keep more info?
    nrows = 0
    for (peptide,iedb_allele,qual_measure) in cursor:
        nrows += 1
        if any(aa not in AAs for aa in peptide): continue
        internal_allele = internal_allele_name(iedb_allele)
        if peptide not in measurements:
            measurements[peptide] = {internal_allele:[qual_measure]}
        elif internal_allele not in measurements[peptide]:
            measurements[peptide][internal_allele] = [qual_measure]
        else:
            measurements[peptide][internal_allele].append(qual_measure)
    print(nrows,'rows','=>',len(measurements),'peptides')

    cursor.close()
    connection.close()
    return measurements

def get_core_hits(measurements, pred, option):
    # TODO: currently just indicates whether there's evidence for a core binding an allele; nothing relative about strength of evidence
    core_hits = collections.defaultdict(set) # core => [allele]
    
    for peptide in measurements:
        if pred:
            pep_epimap = pred.score_protein(peptide)
        for allele in measurements[peptide]:
            data = measurements[peptide][allele]
            if 'Negative' in data:
                if len(set(data))==1: continue # only negative -- skip
                print('conflicting evidence', peptide, allele, data)
                # TODO: specify how to resolve conflicts; for now, err on side of safety (if somebody called it positive, then it's bad)
                # Note that conflicts could also happen at the level of core 9mers from different peptides, but won't be seen as currently coded....
            if option == 'all':
                for i in range(len(peptide)-9): 
                    core_hits[peptide[i:i+9]].add(allele)
            elif option == 'predicted-good':
                a = pred.alleles.find(allele) # TODO: ugh -- details is just a list in order of alleles
                for core in pep_epimap.peptides:
                    if pep_epimap.peptide_score[core].details[a] is not None: # TODO: Propred does this, but not NetMHC (it also doesn't store 9mers, so fix simultaneously)
                        core_hits[core].add(allele)
            else:
                # TODO: implement
                raise Exception('not implemented')
    print(len(core_hits),'cores')
    return core_hits

def setup_parser():
    """Creates an ArgumentParser and sets up all its arguments.
    returns: ArgumentParser"""
    parser = argparse.ArgumentParser(description='Establishes a database of precomputed epitope scores, to enable efficient lookup at design time', 
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog="""Additional notes:
    """)

    # where to get the data
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--iedb_csv', help='name of IEDB-format csv export')
    source.add_argument('--iedb_mysql', help='name of mysql database holding IEDB-format tables')
    # from peptide to consistuent epitope(s)
    parser.add_argument('--cores', help='which 9mer cores of a peptide to include (default: %(default)s)', choices=['all','predicted_good','predicted_best'], default='all')
    
    # epitope predictor
    pred = parser.add_mutually_exclusive_group()
    #pred.add_argument('--pred_db', help='name of database from which to load epitope predictions') # TODO: different name from score.py since output is also db -- confusing?
    #pred.add_argument('--matrix', help='generic epitope predictor matrix filename')
    #pred.add_argument('--netmhcii', action='store_true', help='use netmhcII executable')
    pred.add_argument('--propred', action='store_true', help='use propred matrices')
    # epitope predictor alleles
    #alleles = parser.add_mutually_exclusive_group()
    #alleles.add_argument('--allele_set', help='name of predefined set of alleles', choices=['test', 'greenbaum11', 'paul15', 'all'])
    #alleles.add_argument('--alleles', help='comma-separated list of allele names')
    # epitope predictor parameters
    parser.add_argument('--epi_thresh', help='epitope predictor threshold (default: %(default).2f)', type=int, default=5)
    parser.add_argument('--noncanon', help='how to treat letters other than the 20 canonical AAs (default: %(default)s)', choices=['error', 'silent', 'warn'])
    parser.add_argument('--netmhcii_score', help='type of score to compute (default %(default)s)', choices=['rank','absolute'], default='rank')
    # the database
    parser.add_argument('db', help='name of sqlite3 database to store epitope information (create or augment)')

    return parser

def main(args):
    """Sets up a data-based epitope database based on the parsed args. (Say that 5 times fast.)
    args: ArgumentParser"""
    
    # epitope predictor
    pred = None # default is not to use prediction
    if args.propred:
        pred = Propred.load()
#    elif args.matrix is not None:
#        pred = EpitopePredictorMatrix.load(args.matrix)
#    elif args.netmhcii:
#        pred = NetMHCII(score_type=args.netmhcii_score[0])
#        raise Exception('not implemented yet -- need to extract cores')
#    elif args.db is not None:
#        pred = EpitopeDatabase.for_reading(args.pred_db)
    if pred is not None:
        pred.thresh = args.epi_thresh
        if pred.peptide_length != 9:
            raise Exception('need to predict 9mer cores')

    # alleles
    # TODO: hardcoding for now, due to naming consistency crap
    iedb_alleles = ['HLA-DRB1*01:01'] # TODO: IEDB naming vs. NetMHC naming
    internal_alleles = [internal_allele_name(a) for a in iedb_alleles]
    # TODO: check this....
    iedb_alleles = iedb_alleles + ['HLA-DRA*01:01/'+a[4:] for a in iedb_alleles if a.startswith('HLA-DR')]
    if pred:
        pred.filter_alleles(internal_alleles)
#    if args.allele_set is not None:
#        if args.allele_set not in pred.allele_sets: 
#            raise Exception('allele_set '+args.allele_set+' not supported')
#        pred.filter_alleles(pred.allele_sets[args.allele_set])
#    elif args.alleles is not None:
#        pred.filter_alleles(args.alleles.split(','))

    # load data
    name = 'data from '
    if args.iedb_csv is not None:
        measurements = load_iedb_csv(args.iedb_csv, iedb_alleles)
        name += ' iedb_csv '+args.iedb_csv
    elif args.iedb_mysql is not None:
        measurements = load_iedb_mysql(args.iedb_mysql, iedb_alleles)
        name += ' iedb_mysql '+args.iedb_mysql

    # TODO: when extending to NetMHC, at this point should batch predict for the identified peptides
    # will then need to extract predicted core
    
    # which cores to store        
    if args.cores == 'all':
        name += ' all cores'
    elif args.cores == 'predicted_good' or args.cores == 'predicted_best':
        if pred is None:
            raise Exception('if going to use predicted cores (--cores predicted_*), must specify an epitope predictor')
        # TODO: NetMHC only returns one predicted core, so disallow good or at least print warning
        name += ' '+args.cores+ ' '+pred.name
    core_hits = get_core_hits(measurements, pred, args.cores)
    
    # make scores for the cores
    cores = list(core_hits.keys())
    # TODO: current binary score (is there a hit or not); could scale based on amount of evidence, # alleles, whatever
    scores = [EpitopeScore(1, [1 if a in hits else None for a in internal_alleles]) for hits in core_hits.values()]
    epimap = EpitopeMap(9, internal_alleles, cores, scores)
        
    db = EpitopeDatabase.for_writing(args.db, name, internal_alleles, 9)
    db.save_scores(epimap)
        
if __name__ == '__main__':
    parser = setup_parser()
    args = parser.parse_args()
    main(args)
    # TODO: protected main
#    try:
#        main(args)
#    except Exception as e:
#        print('ERROR', e)
#        print()
#        parser.print_usage()
