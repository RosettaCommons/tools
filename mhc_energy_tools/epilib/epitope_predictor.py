#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Base classes for epitope prediction

@author: Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu 
"""

class EpitopeScore (object):
    """A score for a peptide, both an overall value and details for predicted alleles."""
    def __init__(self, value=0, details=[]):
        self.value = value
        self.details = details

class EpitopeMap (object):
    """Scores for peptides."""
    def __init__(self, peptide_length, alleles, peptides, scores):
        self.peptide_length = peptide_length
        self.alleles = alleles
        self.peptides = peptides
        self.scores = scores
        self.p2s = dict(zip(peptides, scores))
    
    def peptide_score(self, peptide):
        """Returns the score for the peptide"""
        return self.p2s[peptide]
    
    def total_score(self):
        """Sum of scores over the peptides"""
        return sum(score.value for score in self.scores)

class EpitopePredictor (object):
    """Provides the ability to predict epitope scores for peptides (or a whole protein) against MHC alleles."""
    
    # predefined sets of alleles that can be accessed by name for this predictor
    # { String set_name : [ String allele_name ] }
    allele_sets = {}
    
    def __init__(self, name, alleles=[], peptide_length=9):
        self.name = name
        self.alleles = alleles
        self.peptide_length = peptide_length
 
    def filter_alleles(self, selected_alleles):
        """Subsets the predictor's alleles to those listed, making sure they're supported."""
        missing = set(selected_alleles) - set(self.alleles)
        if len(missing)>0:
            raise Exception('no such allele(s) '+','.join(missing))
        self.alleles = selected_alleles
        
    def score_peptide(self, pep):
        """EpitopeScore for a single peptide (as a String)."""
        raise Exception('score_peptide needs to be implemented by subclass')
    
    def score_peptides(self, peptides):
        """EpitopeScores for a list of peptides (each a String). May be more efficient in batch than one-by-one."""
        return EpitopeMap(self.peptide_length, self.alleles,
                          peptides, [self.score_peptide(pep) for pep in peptides])
    
    def score_protein(self, seq):
        """EpitopeMap for a whole protein (as a String)."""
        peptides = [seq[i:i+self.peptide_length] for i in range(len(seq)-self.peptide_length+1)]
        return EpitopeMap(self.peptide_length, self.alleles,
                          peptides, [self.score_peptide(pep) for pep in peptides])