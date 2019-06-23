#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parsing experimental epitope data, currently with support for IEDB-downloaded files or the complete IEDB database

@author: Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu 
"""

"""Notes during development
* example query https://help.iedb.org/hc/en-us/articles/114094146451-Select-all-MHC-binding-assays-for-class-II-human-epitopes-with-a-SQL-queryAP

select m.mhc_bind_id into outfile "~/Downloads/bind.csv" from mhc_bind m where m.mhc_allele_name in ('HLA-DRB1*01:01', 'HLA-DRA*01:01/DRB1*01:01') order by m.mhc_bind_id;
select m.mhc_elution_id into outfile "~/Downloads/elute.csv" from mhc_elution m where m.mhc_allele_name in ('HLA-DRB1*01:01','HLA-DRA*01:01/DRB1*01:01') order by m.mhc_elution_id;
-> combining these gives almost the same as what I get via web query (just restricting HLA type to 0101); web query has two additional sets of epitopes, apparently from two additional publications (recently added?)

* fields http://curationwiki.iedb.org/wiki/index.php/Data_Field_Descriptions
"""

import collections, csv, os, re

from epilib.sequence import AAs

class ExptData (object):
    """Manages experimental epitope data.
    Established as a base class to support future expansion possibilities (e.g., simple csv format)."""
 
    # predefined sets of alleles that can be accessed by name for this predictor
    # { String set_name : [ String allele_name ] }
    # TODO: copied from EpitopePredictor -- generalize somehow?
    allele_sets = {}

    def __init__(self, name, std_alleles, measurements, scores):
        """name: arbitrary string
        std_alleles: allele name => standardized allele name
        measurements: peptide => allele => [data], where data is implementation-dependent, but somehow measures the allele:peptide binding
        scores: peptide => allele => number, compiled from measurements into a score"""
        self.name = name
        self.std_alleles = std_alleles
        self.measurements = measurements
        self.scores = scores
    

class IEDBData (ExptData):
    """Manages experimental epitope data downloaded from the IEDB."""
      
    # names as IEDB wants them
    allele_sets = {
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
        
        # Southwood S, Sidney J, Kondo A, del Guercio MF, Appella E, Hoffman S, Kubo RT, Chesnut RW, Grey HM, Sette A.
        # Several common HLA-DR types share largely overlapping peptide binding repertoires
        # J Immunol. 1998 Apr;160:3363-73
        # https://www.ncbi.nlm.nih.gov/pubmed/9531296
        'southwood98':
            ['HLA-DRB1*01:01','HLA-DRB1*03:01','HLA-DRB1*04:01','HLA-DRB1*07:01','HLA-DRB1*08:01','HLA-DRB1*11:01','HLA-DRB1*13:01','HLA-DRB1*15:01']
        }

    # TODO: copied and modified from NetMHC -- generalize?
    # in particular, use the internal "std_name" to reconcile naming    
    iedb2std = {} # cache for std_name
    @staticmethod
    def std_name(iedb_name):
        """Converts IEDB name into something suitable for a column name in the database.
        - drops "HLA-"
        - drops ':'
        - converts other non-alphanumeric characters to _
        Caches the conversion."""
        if iedb_name in IEDBData.iedb2std: return IEDBData.iedb2std[iedb_name]
        name = iedb_name \
            .replace('HLA-','') \
            .replace(':', '')
        name = re.sub('[^A-Za-z0-9_]', '_', name)
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
    def from_csv(filename, alleles=None, allele_set=None, allele_file=None):
        """Loads from IEDB-formatted download file.
        If given, restricts to just the specified alleles (by list of names, named set, or file); else imports all.
        Special allele_set name 'hlaII' means from file 'data/iedb-hlaII.csv'."""

        restr_synonyms = None # restricted query allele name => synonym for merging
        if allele_set == 'hlaII': restr_synonyms = IEDBData.load_allele_file(os.path.dirname(__file__) + '/../data/iedb-hlaII.csv')
        elif allele_set is not None: restr_synonyms = IEDBData.loaded_synonyms(IEDBData.allele_sets[allele_set])
        elif alleles is not None: restr_synonyms = IEDBData.loaded_synonyms(alleles)
        elif allele_file is not None: restr_synonyms = IEDBData.load_allele_file(allele_file)
        std_alleles = dict((a,IEDBData.std_name(s)) for (a,s) in restr_synonyms.items())

        measurements = collections.defaultdict(lambda:collections.defaultdict(list)) # peptide => allele => [qual_measure]
        std_alleles = {}
        with open(filename) as infile:
            # IEDB splits the header over two lines, so concatenate them as use them as the field names
            head1 = infile.readline().strip().split(',')
            head2 = infile.readline().strip().split(',')
            nrows = 0
            for row in csv.DictReader(infile, fieldnames=[h1+' '+h2 for h1,h2 in zip(head1,head2)]):
                nrows += 1
                peptide = row['Epitope Description']
                allele = row['MHC Allele Name']
                if restr_synonyms is not None and allele not in restr_synonyms:
                    #print('skipping allele',allele)
                    continue
                # only deal with peptides containing at least 9 amino acids and nothing funky
                if len(peptide)<9 or any(aa not in AAs for aa in peptide):
                    print('skipping',peptide)
                    continue
                if allele not in std_alleles:
                    std_alleles[allele] = IEDBData.std_name(allele)
                qual_measure = row['Assay Qualitative Measure']
                measurements[peptide][std_alleles[allele]].append(qual_measure)
        print(nrows,'rows','=>',len(measurements),'peptides')
        return IEDBData('iedb csv '+filename, std_alleles, measurements, IEDBData.aggregate_measurements(measurements))

    @staticmethod
    def loaded_synonyms(alleles):
        synonyms = IEDBData.load_allele_file(os.path.dirname(__file__) + '/../data/iedb-hlaII.csv')
        return dict((a, synonyms[a] if a in synonyms else a) for a in alleles)

    @staticmethod
    def download_mysql(dbname, user='root', pw=None):
        """Download IEDB database into local mysql database (blowing away whatever's there now)."""
        mysql = 'mysql -u '+user
        if pw is not None: mysql += ' -p'+pw
        # http://www.iedb.org/database_export_v3.php
        errcode = os.system(mysql + ' -e "drop database if exists '+dbname+'; create database '+dbname+';"')
        if (errcode):
            raise Exception('Error re-creating empty database ' + dbname + ', with error code ' + str(errcode))
        errcode = os.system('curl http://www.iedb.org/downloader.php?file_name=doc/iedb_public.sql.gz | gunzip -c | '+ mysql + ' iedb ')
        if (errcode):
            raise Exception('Error downloading IEDB or populating local database ' + dbname + ', with error code ' + str(errcode))

    @staticmethod
    def from_mysql(dbname, assay_binding_filter, assay_elution_filter, alleles=None, allele_set=None, allele_file=None, user='root', pw=None):
        """Loads from local mysql database downloaded from IEDB.
        Either or both of data from the mhc_bind table and/or the mhc_elution table.
        A list of alleles, an allele_set name, or an allele_file must be given.
        Special allele_set name 'hlaII' means everything that is returned in a query for human class II."""

        import mysql.connector
        try:
            connection = mysql.connector.connect(database=dbname, user=user, password=pw) # TODO: server, other options? https://dev.mysql.com/doc/connector-python/en/connector-python-connectargs.html
            cursor = connection.cursor()
        except mysql.connector.Error as err:
            print('Error accessing local mysql database ' + dbname)
            raise(err)

        synonyms = None # query allele name => synonym for merging
        if allele_set == 'hlaII': synonyms = IEDBData.query_hlaII(cursor)
        elif allele_set is not None: synonyms = IEDBData.query_synonyms(cursor, IEDBData.allele_sets[allele_set])      
        elif alleles is not None: synonyms = IEDBData.query_synonyms(cursor, alleles)
        elif allele_file is not None: synonyms = IEDBData.load_allele_file(allele_file)
        else: raise Exception('please specify alleles using one of the available mechanisms')
        std_alleles = dict((a,IEDBData.std_name(s)) for (a,s) in synonyms.items())
        
        measurements = collections.defaultdict(lambda:collections.defaultdict(list)) # peptide => allele => [qual_measure]
        if assay_binding_filter != 'none':
            print('loading mhc ligand binding data')
            IEDBData.load_table(cursor, measurements, 'mhc_bind', assay_binding_filter, std_alleles)
        if assay_elution_filter != 'none':
            print('loading mhc ligand elution data')
            IEDBData.load_table(cursor, measurements, 'mhc_elution', assay_elution_filter, std_alleles)
        # TODO: t cell data?
        print('=>',len(measurements),'peptides')

        cursor.close()
        connection.close()

        return IEDBData('iedb mysql '+dbname, std_alleles, measurements, IEDBData.aggregate_measurements(measurements))

    @staticmethod
    def query_synonyms(cursor, alleles):
        """Uses the cursor to find synonyms for the given alleles.
        Currently for DRB-only alleles (no mutations) also finds those with DRA*01:01 (no mutations) specified (this is what the IEDB web site does).
        Returns allele => synonym (where synonym is allele itself for most)"""

        synonyms = dict((a,a) for a in alleles)
        drbs = ['"%s"' % (a,) for a in alleles if a.startswith('HLA-DRB')]
        print(drbs)
        if not drbs: return alleles # If there are no drbs, we don't need to worry about alleles.  Just give back the input alleles.
        cursor.execute('select displayed_restriction,chain_ii_name from mhc_allele_restriction where chain_i_name="HLA-DRA*01:01" and chain_i_mutation is NULL and chain_ii_name in (%s) and chain_ii_mutation is NULL' % (','.join(drbs), ))
        # TODO: error handling

        for (allele,synonym) in cursor:
            # The mysql connector sometimes (but not always) returns bytearrays instead of strings.  Decode them to strings only if necessary.
            if type(allele) == bytearray: allele = allele.decode()
            if type(synonym) == bytearray: synonym = synonym.decode()
            synonyms[allele] = synonym

        for a_s in sorted(synonyms.items()): print(a_s)
        
        return synonyms
        
    @staticmethod
    def query_hlaII(cursor):
        """Uses the cursor to query HLA-II alleles. 
        Identifies synonyms for DRA/DRB as in get_synonyms.
        Returns allele => synonym (where synonym is allele itself for most)"""

        cursor.execute('select displayed_restriction,chain_i_name,chain_i_mutation,chain_ii_name,chain_ii_mutation from mhc_allele_restriction where class="II" and organism_ncbi_tax_id="9606" order by displayed_restriction')
        # TODO: error handling

        synonyms = {}     
        for (allele,a_chain,a_mut,b_chain,b_mut) in cursor:
            # The mysql connector sometimes (but not always) returns bytearrays instead of strings.  Decode them to strings only if necessary.
            if type(allele) == bytearray: allele = allele.decode()
            if type(a_chain) == bytearray: a_chain = a_chain.decode()
            if type(b_chain) == bytearray: b_chain = b_chain.decode()
            synonyms[allele] = allele
            if a_chain=='HLA-DRA*01:01' and a_mut is None and b_chain.startswith('HLA-DRB') and b_mut is None:
                synonyms[allele] = b_chain
                
        return synonyms

    @staticmethod
    def load_allele_file(filename):
        """Loads the allele names (and optional synonyms) from the csv file.
        Format: header row, then one row per allele name (first column) with optional synonym (second column) and epitope predictor allele name (third column, currently ignored)
        Returns allele => synonym (itself if unspecified)"""
        
        synonyms = {}
        with open(filename,'r') as infile:
            infile.readline() # hedaer
            for row in csv.reader(infile):
                if row[1] is not None: synonyms[row[0]] = row[1]
                else: synonyms[row[0]] = row[0]
                # TODO: get epipred name
                
        return synonyms
                
    @staticmethod
    def load_table(cursor, measurements, table, assay_filter, std_names):
        # TODO: allow filtering on specific assay details (currently just ignoring assay_filter)

        # Note: while the IEDB example suggests to get the linear_peptide_seq, that doesn't include modifications, whereas description seems to
        # (e.g., id 95191 is 'AADAAAKAAAAAAA + MCM(3)' and a comment says that 'The third residue of the epitope is benzylaspartic acid')
        # TODO: filter out modified peptides as part of the query?
        query = """
        select e.description as peptide, m.mhc_allele_name as allele, m.as_char_value as qual_measure 
        from %s m, curated_epitope ce, epitope_object eo, epitope e
        where m.curated_epitope_id=ce.curated_epitope_id and ce.e_object_id = eo.object_id and eo.epitope_id=e.epitope_id
        and length(e.linear_peptide_seq)>=9
        and m.mhc_allele_name in (%s)""" % (table, ','.join('"'+a+'"' for a in std_names))
        # print(query)
        
        cursor.execute(query)
        # TODO: error handling
        
        nrows = 0
        for (peptide,iedb_allele,qual_measure) in cursor:
            nrows += 1
            
            # The mysql connector sometimes (but not always) returns bytearrays instead of strings.  Decode them to strings only if necessary.
            if type(peptide) == bytearray: peptide = peptide.decode()
            if type(iedb_allele) == bytearray: iedb_allele = iedb_allele.decode()
            if type(qual_measure) == bytearray: qual_measure = qual_measure.decode()
            
            # only deal with peptides containing nothing funky
            # (in contrast to csv, the length condition was handled during the query; TODO: if it's of interest to print it out, modify accordingly)
            if any(aa not in AAs for aa in peptide):
                print('skipping',peptide)
                continue
            measurements[peptide][std_names[iedb_allele]].append(qual_measure)
        print(nrows,'rows')

    @staticmethod
    def export_hlaII(dbname, filename, user='root', pw=None):
        """Internal method (not currently available from command-line, but could be) to extract HLA-II allele names into a csv file."""

        import mysql.connector
        try:
            connection = mysql.connector.connect(database=dbname, user=user, password=pw) # TODO: server, other options? https://dev.mysql.com/doc/connector-python/en/connector-python-connectargs.html
            cursor = connection.cursor()
        except mysql.connector.Error as err:
            print('Error accessing local mysql database ' + dbname)
            raise(err)

        synonyms = IEDBData.query_hlaII(cursor)
        
        with open(filename,'w') as outfile:
            outcsv = csv.writer(outfile)
            outcsv.writerow(['allele','synonym'])
            
            for allele_synonym in sorted(synonyms.items()):
                outcsv.writerow(allele_synonym)

        cursor.close()
        connection.close()
