#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parsing experimental epitope data, currently with support for IEDB-downloaded files or the complete IEDB database

@author: Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu 
"""

import collections, csv, os

from epilib.sequence import AAs

class ExptData (object):
    """Manages experimental epitope data.
    Established as a base class to support future expansion possibilities (e.g., simple csv format)."""
 
    # predefined sets of alleles that can be accessed by name for this predictor
    # { String set_name : [ String allele_name ] }
    # TODO: copied from EpitopePredictor -- generalize somehow?
    allele_sets = {}

    @staticmethod
    def std_name(name):
        """Converts name as given in the file to a 'standard' name used in the rest of the library"""
        return name

    def __init__(self, name, std_alleles, measurements, scores):
        """name: arbitrary string
        std_alleles: list of std allele names for which there are measurements
        measurements: peptide => allele => [data], where data is implementation-dependent, but somehow measures the allele:peptide binding
        scores: peptide => allele => number, compiled from measurements into a score"""
        self.name = name
        self.std_alleles = std_alleles
        self.measurements = measurements
        self.scores = scores
    

class IEDBData (ExptData):
    """Manages experimental epitope data downloaded from the IEDB."""
      
    # names as IEDB wants them
    # also adding in DRA partner allele since some data includes it
    # TODO: check legit
    def add_dra(allele_sets):
        return dict((name, alleles+['HLA-DRA*01:01/'+a[4:] for a in alleles if a.startswith('HLA-DR')]) for (name,alleles) in allele_sets.items())
    allele_sets = add_dra({
        'test': 
            ['HLA-DRB1*01:01'],
            
        # TODO: check these
        
        # Greenbaum J, Sidney J, Chung J, Brander C, Peters B, Sette A.
        # Functional classification of class II human leukocyte antigen (HLA) molecules reveals seven different supertypes and a surprising degree of repertoire sharing across supertypes
        # Immunogenetics. 2011 Jun;63(6):325-35
        # https://www.ncbi.nlm.nih.gov/pubmed/21305276
        'greenbaum11':
            ['HLA-DRB1*01:01', 'HLA-DRB1*03:01', 'HLA-DRB1*04:01', 'HLA-DRB1*04:05', 'HLA-DRB1*07:01', 'HLA-DRB1*08:02', 'HLA-DRB1*09:01', 'HLA-DRB1*11:01', 'HLA-DRB1*12:01', 'HLA-DRB1*13:02', 'HLA-DRB1*15:01',
             'HLA-DRB3*01:01', 'HLA-DRB3*02:02', 'HLA-DRB4*01:01', 'HLA-DRB5*01:01',
             'HLA-DPA1*01:03/DPB1*02:01', 'HLA-DPA1*01:03/DPB1*04:01', 'HLA-DPA1*02:01/DPB1*01:01', 'HLA-DPA1*02:01/DPB1*05:01', 'HLA-DPA1*02:01/DPB1*14:01', 'HLA-DPA1*03:01/DPB1*04:02',
             'HLA-DQA1*01:01/DQB1*05:01', 'HLA-DQA1*01:02/DQB1*06:02', 'HLA-DQA1*03:01/DQB1*03:02', 'HLA-DQA1*04:01/DQB1*04:02', 'HLA-DQA1*05:01/DQB1*02:01', 'HLA-DQA1*05:01/DQB1*03:01'],

        # Paul S, Lindestam Arlehamn CS, Scriba TJ, Dillon MB, Oseroff C, Hinz D, McKinney DM, Carrasco Pro S, Sidney J, Peters B, Sette A.
        # Development and validation of a broad scheme for prediction of HLA class II restricted T cell epitopes
        # J Immunol Methods. 2015 Jul;422:28-34
        # https://www.ncbi.nlm.nih.gov/pubmed/25862607
        'paul15':
            ['HLA-DRB1*03:01', 'HLA-DRB1*07:01', 'HLA-DRB1*15:01', 'HLA-DRB3*01:01', 'HLA-DRB3*02:02', 'HLA-DRB4*01:01', 'HLA-DRB5*01:01'],
        })

    # TODO: copied and modified from NetMHC -- generalize?
    # in particular, use the internal "std_name" to reconcile naming    
    iedb2std = {} # cache for std_name
    @staticmethod
    def std_name(iedb_name):
        """Converts IEDB name into something suitable for a column name in the database.
        - drops "HLA-"
        - converts "-", '*', and '/' to "_"
        - drops ':'
        - gets rid of DRA paired with DRB [for uniformity, collapsing; see TODO below]
        Caches the conversion."""
        if iedb_name in IEDBData.iedb2std: return IEDBData.iedb2std[iedb_name]
        name = iedb_name \
            .replace('HLA-','') \
            .translate(str.maketrans({'/':'_', '-':'_', '*':'_', ':':None})) \
            .replace('DRA_0101_','') 
            # TODO: I think we want to collapse the DRA down, but double check or allow option?
        IEDBData.iedb2std[iedb_name] = name
        return name

    @staticmethod
    def aggregate_measurements(measurements):
        """From measurements to scores, for constructor"""
        
        # TODO: specify how to resolve conflicts; for now, err on side of safety: if somebody called it positive, then it's bad
        # note that conflicts could also happen at the level of core 9mers from different peptides, but won't be seen as currently coded....

        scores = collections.defaultdict(dict) # peptide => allele => score
        
        for peptide in measurements:
            for allele in measurements[peptide]:
                data = measurements[peptide][allele]
                if 'Negative' in data:
                    if len(set(data))==1: continue # only negative -- skip
                    print('conflicting evidence', peptide, allele, data)
                scores[peptide][allele] = 1
                
        return scores

    @staticmethod
    def from_csv(filename, alleles=None, allele_set=None):
        """Loads from IEDB-formatted download file.
        If given, restricts to just the specified alleles or named allele_set; else imports all."""
        if allele_set is not None: alleles = IEDBData.allele_sets[allele_set]
        measurements = collections.defaultdict(lambda:collections.defaultdict(list)) # peptide => allele => [qual_measure]
        std_alleles = set()
        with open(filename) as infile:
            # IEDB splits the header over two lines, so concatenate them as use them as the field names
            head1 = infile.readline().strip().split(',')
            head2 = infile.readline().strip().split(',')
            nrows = 0
            for row in csv.DictReader(infile, fieldnames=[h1+' '+h2 for h1,h2 in zip(head1,head2)]):
                nrows += 1
                peptide = row['Epitope Description']
                allele = row['MHC Allele Name']
                if alleles is not None and allele not in alleles:
                    #print('skipping allele',allele)
                    continue
                # only deal with peptides containing at least 9 amino acids
                if len(peptide)<9 or any(aa not in AAs for aa in peptide): continue
                std_allele = IEDBData.std_name(allele)
                std_alleles.add(std_allele)
                qual_measure = row['Assay Qualitative Measure']
                measurements[peptide][std_allele].append(qual_measure)
        print(nrows,'rows','=>',len(measurements),'peptides')
        return IEDBData('iedb csv '+filename, std_alleles, measurements, IEDBData.aggregate_measurements(measurements))

    @staticmethod
    def download_mysql(dbname, user='root', pw=None):
        """Download IEDB database into local mysql database (blowing away whatever's there now)."""
        mysql = 'mysql -u '+user
        if pw is not None: mysql += ' -p'+pw
        # http://www.iedb.org/database_export_v3.php
        os.system(mysql + ' -e "drop database if exists '+dbname+'; create database '+dbname+';"')
        os.system('curl http://www.iedb.org/downloader.php?file_name=doc/iedb_public.sql.gz | gunzip -c | '+ mysql + ' iedb ')

    @staticmethod
    def from_mysql(dbname, alleles=None, allele_set=None, user='root', pw=None):
        """Loads from local mysql database downloaded from IEDB.
        A list of alleles or an allele_set name must be given."""
        if allele_set is not None: alleles = IEDBData.allele_sets[allele_set]       
        if alleles is None: raise Exception('please specify alleles or allele_set')
        std_alleles = set(IEDBData.std_name(a) for a in alleles) # duplicates due to DRA
        
        """Notes during development
        * example query https://help.iedb.org/hc/en-us/articles/114094146451-Select-all-MHC-binding-assays-for-class-II-human-epitopes-with-a-SQL-queryAP

        select m.mhc_bind_id into outfile "~/Downloads/bind.csv" from mhc_bind m where m.mhc_allele_name in ('HLA-DRB1*01:01', 'HLA-DRA*01:01/DRB1*01:01') order by m.mhc_bind_id;
        select m.mhc_elution_id into outfile "~/Downloads/elute.csv" from mhc_elution m where m.mhc_allele_name in ('HLA-DRB1*01:01','HLA-DRA*01:01/DRB1*01:01') order by m.mhc_elution_id;
        -> combining these gives almost the same as what I get via web query (just restricting HLA type to 0101); web query has two additional sets of epitopes, apparently from two additional publications (recently added?)

        * fields http://curationwiki.iedb.org/wiki/index.php/Data_Field_Descriptions
        """

        # TODO: currently just binding data; also get elution from separate table (as in note above)
        # TODO: t cell data?
        # TODO: allow filtering on assay details
        query = ('select e.linear_peptide_seq as peptide, m.mhc_allele_name as allele, m.as_char_value as qual_measure '
                 'from mhc_bind m, curated_epitope ce, epitope_object eo, epitope e '
                 'where m.curated_epitope_id=ce.curated_epitope_id and ce.e_object_id = eo.object_id and eo.epitope_id=e.epitope_id '
                 'and length(e.linear_peptide_seq)>=9 '
                 'and m.mhc_allele_name in ') + '(' + ','.join('"'+a+'"' for a in alleles) + ')'
        print(query)
        
        import mysql.connector
        connection = mysql.connector.connect(database=dbname, user=user, password=pw) # TODO: server, ... ? https://dev.mysql.com/doc/connector-python/en/connector-python-connectargs.html
        cursor = connection.cursor()
        cursor.execute(query)
        # TODO: error handling
        
        measurements = collections.defaultdict(lambda:collections.defaultdict(list)) # peptide => allele => [qual_measure]
        nrows = 0
        for (peptide,iedb_allele,qual_measure) in cursor:
            nrows += 1
            if any(aa not in AAs for aa in peptide): continue
            measurements[peptide][IEDBData.std_name(iedb_allele)].append(qual_measure)
        print(nrows,'rows','=>',len(measurements),'peptides')
    
        cursor.close()
        connection.close()

        return IEDBData('iedb mysql '+dbname, std_alleles, measurements, IEDBData.aggregate_measurements(measurements))
