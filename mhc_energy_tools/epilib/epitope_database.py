#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Maintenance of a database of epitope scores, caching results from an epitope predictor for efficient lookup in design.

@author: Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu 
"""

import sqlite3, os.path
from epilib.epitope_predictor import EpitopePredictor, EpitopeScore, EpitopeMap

class EpitopeDatabase (EpitopePredictor):  
    """Manages the database and serves itself as an EpitopePredictor"""
    
    def __init__(self, conn, cursor, name, alleles=[], peptide_length=9, is_new=False):
        super().__init__(name, alleles, peptide_length)
        self.conn = conn
        self.cursor = cursor
        self.is_new = is_new  # creating, rather than adding to existing
    
    @classmethod
    def for_reading(cls, filename, handle_unseen='w', unseen_score=100):
        """Opens an existing database for reading (i.e., acting as an EpitopePredictor).
        Handle_unseen indicates whether to warn ('w'), return the unseen_score ('s') or else raise an exception."""
        if not os.path.exists(filename): raise Exception("database '" + filename+"' doesn't exist")
        conn = sqlite3.connect(filename)
        conn.row_factory = sqlite3.Row # so can access columns by name
        cursor = conn.cursor()
        # get predictor info from db
        try:
            alleles = [row['name'] for row in cursor.execute('select * from alleles')]
            row = cursor.execute('select value from meta where name="peptide_length"').fetchone()
            peptide_length = int(row['value'])
        except:
            raise Exception("database '" + filename+"' not properly constructed")
        db = cls(conn, cursor, filename, alleles, peptide_length)
        db.handle_unseen = handle_unseen
        db.unseen_score = unseen_score
        return db
    
    @classmethod
    def for_writing(cls, filename, source, alleles=[], peptide_length=9):
        """Opens a database for storing epitope data.
        If the database already exists, the alleles and peptide_length must be the same."""
        conn = sqlite3.connect(filename)
        conn.row_factory = sqlite3.Row # so can access columns by name
        cursor = conn.cursor()
        row = cursor.execute('select count(*) as value from sqlite_master where type="table" and name="meta" or name="alleles" or name="epitopes"').fetchone();
        count = int(row['value'])
        if count==0:
            # new database; set up
            cursor.execute('create table meta (name text primary key, value text)')
            cursor.execute('insert into meta values (?, ?)', ('source', source))
            cursor.execute('insert into meta values (?, ?)', ('peptide_length', peptide_length))
            cursor.execute('create table alleles (name text primary key)')
            cursor.executemany('insert into alleles values (?)', [(a,) for a in alleles])
            cursor.execute('create table epitopes (peptide text primary key, score real, ' + ','.join(alleles) + ')')
            conn.commit()
        elif count==3:
            # existing database; make sure parameters match
            errors = []
            row = cursor.execute('select value from meta where name="source"').fetchone()
            # seems like this might not always be appropriate, e.g., when combining different sources of data
            # so now 'source' is essentially just for informational purposes
#            if row['value'] != source:
#                errors.append('mismatched source (%s vs. %s)' % (row['value'], source))
            row = cursor.execute('select value from meta where name="peptide_length"').fetchone()
            if int(row['value']) != peptide_length: 
                errors.append('mismatched peptide length (%d vs. %d)' % (int(row['value']), peptide_length))
            db_alleles = set(row['name'] for row in cursor.execute('select * from alleles'))
            if db_alleles != set(alleles):
                errors.append('mismatched alleles (%d vs. %d of them)' % (len(db_alleles), len(set(alleles))))
            if len(errors)>0:
                raise Exception('database exists but established with different parameters: ' + '; '.join(errors))
        else:
            raise Exception('database exists but schema is wrong')

        db = cls(conn, cursor, filename, alleles, peptide_length, is_new=count==0)

        # note: assuming the predictor will return the same value every time, so just ignore
        db.save_stmt = 'insert or ignore into epitopes values (?,?,' + ','.join('?' for a in alleles) + ')'

        return db

    def has_peptide(self, peptide):
        """Is the peptide in the db at all?"""
        try:
            row = self.cursor.execute('select score from epitopes where peptide=?', (peptide,)).fetchone()
            return row is not None
        except KeyError:
            return False
             
    def score_peptide(self, peptide):
        try:
            row = self.cursor.execute('select * from epitopes where peptide=?', (peptide,)).fetchone()
            if row is not None: return EpitopeScore(row['score'], [row[a] for a in self.alleles])
        except KeyError as e:
            row = None
        if self.handle_unseen == 'w':
            print('warning: unscored peptide '+peptide)
            return EpitopeScore(self.unseen_score, [])
        elif self.handle_unseen == 's':
            return EpitopeScore(self.unseen_score, [])
        else:
            raise Exception('unscored peptide '+peptide)
                
    def save_scores(self, epimap):
        """Saves the scores from the EpitopeMap to the database"""
        # epimap.report()
        self.cursor.executemany(self.save_stmt, [[p, epimap.peptide_score(p).value] + epimap.peptide_score(p).details for p in epimap.peptides])
        self.conn.commit()

    def load_scores(self):
        """Returns an EpitopeMap with all the scores"""
        scores = {}
        for row in self.cursor.execute('select * from epitopes'):
            scores[row['peptide']] = EpitopeScore(row['score'], [row[a] for a in self.alleles])
        return EpitopeMap(self.peptide_length, self.alleles, scores.keys(), scores.values())
        