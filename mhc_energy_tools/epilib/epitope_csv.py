#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Maintenance of epitope scores, caching results from an epitope predictor for efficient lookup in design.
Essentially a stripped-down version of epitope_database, stored in a csv-format file and read into memory.
Didn't bother coming up with a baseclass to rule them both, though.
The scores are one per row, with the first two columns giving the peptide and overall score, and the remaining ones the details for the alleles.

@author: Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu 
"""

import csv, os.path
from epilib.epitope_predictor import EpitopePredictor, EpitopeScore, EpitopeMap

class EpitopeCSV (EpitopePredictor):  
    """Manages the csv and serves itself as an EpitopePredictor"""
    
    def __init__(self, name, scores=None, alleles=[], peptide_length=9):
        super().__init__(name, alleles, peptide_length)
        if scores is None: self.scores = {}
        else: self.scores = scores
    
    @classmethod
    def for_reading(cls, filename, handle_unseen='w', unseen_score=100):
        """Opens an existing csv for reading (i.e., acting as an EpitopePredictor).
        Reads metadata before 'start', and then reads all the scores into a dictionary for lookup.
        Handle_unseen indicates whether to warn ('w'), return the unseen_score ('s') or else raise an exception."""
        if not os.path.exists(filename): raise Exception(filename+"' doesn't exist")
        with open(filename, 'r') as f:
            # header row -- rigid format
            header = f.readline().strip().split(',')
            if len(header) < 2 or header[0] != 'peptide' or header[1] != 'score':
                raise Exception('bad header row in csv file '+filename)
            alleles = header[2:]
            # get scores
            scores = {}
            for row in csv.reader(f):
                try:
                    peptide = row[0]
                    score = float(row[1])
                    # assumes empty cells are there, not missing
                    details = [0 if v=='' else float(v) for v in row[2:len(header)]]
                    if peptide in scores:
                        # TODO: raise exception?
                        print('duplicate peptide '+peptide+' -- ignoring second one')
                    else:
                        scores[peptide] = EpitopeScore(score, details)
                except:
                    raise Exception('bad score row in csv '+filename+':' + str(row))
        lens = set(len(p) for p in scores)
        if len(lens)==0:
            # TODO: raise exception?
            print('empty csv file '+filename)
            peptide_length = 0
        elif len(lens)==1:
            peptide_length = lens.pop()
        else:
            raise Exception('mismatched peptide lengths not yet handled '+str(lens))
        obj = cls(filename, scores, alleles, peptide_length)
        obj.handle_unseen = handle_unseen
        obj.unseen_score = unseen_score
        return obj
    
    @classmethod
    def for_writing(cls, filename, alleles, peptide_length=9):
        """Opens a csv for storing epitope data."""
        file = open(filename,'w')
        csvf = csv.writer(file, lineterminator=os.linesep)
        # header
        csvf.writerow(['peptide','score']+alleles)
            
        obj = cls(filename, alleles=alleles, peptide_length=peptide_length)
        obj.file = file
        obj.csvf = csvf
        return obj

    def has_peptide(self, peptide):
        """Is the peptide in the db at all?"""
        return peptide in self.scores
             
    def score_peptide(self, peptide):
        if peptide in self.scores:
            return self.scores[peptide]
        elif self.handle_unseen == 'w':
            print('warning: unscored peptide '+peptide)
            return EpitopeScore(self.unseen_score, [])
        elif self.handle_unseen == 's':
            return EpitopeScore(self.unseen_score, [])
        else:
            raise Exception('unscored peptide '+peptide)
                
    def save_scores(self, epimap):
        """Saves the scores from the EpitopeMap to the csv file"""
        for peptide in epimap.peptides:
            score = epimap.peptide_score(peptide)
            self.csvf.writerow([peptide,score.value]+score.details)
            
    def load_scores(self):
        """Returns an EpitopeMap with all the scores"""
        return EpitopeMap(self.peptide_length, self.alleles, self.scores.keys(), self.scores.values())
