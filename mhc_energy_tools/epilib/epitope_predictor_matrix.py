#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Bare bones epitope predictor based on allele-specific position-weight matrices like Propred

To use a matrix in the Rosetta mhc_epitope database, set enviornment variable ROSETTA to the root directory of the Rosetta installation (containing main/database/...)

@author: Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu 
"""

import os
from epilib.epitope_predictor import EpitopePredictor, EpitopeScore

class AlleleMatrix (object):
    def __init__(self, name, threshes, profile):
        self.name = name
        self.threshes = threshes
        self.profile = profile

    def is_hit(self, pep, thresh):
        total = sum(self.profile[i][pep[i]] for i in range(len(pep)))
        return total >= self.threshes[thresh]
    
    def score(self, pep, thresh):
        total = sum(self.profile[i][pep[i]] for i in range(len(pep)))
        for thresh in range(1, thresh+1):
            if total >= self.threshes[thresh]: return thresh
        return None

class EpitopePredictorMatrix (EpitopePredictor):
    def __init__(self, name, matrices=[], peptide_length=9, unknown_aa='w', thresh=5):
        super().__init__(name, alleles=[m.name for m in matrices], peptide_length=peptide_length)
        self.matrices = matrices
        self.unknown_aa = unknown_aa
        self.thresh = thresh
    
    def filter_alleles(self, selected_alleles):
        a2m = dict((m.name, m) for m in self.matrices)
        selected_matrices = []
        missing = []
        for a in selected_alleles:
            if a not in a2m: missing.append(a)
            else: selected_matrices.append(a2m[a])
        if len(missing)>0:
            raise Exception('missing matrix/ces '+','.join(missing))
        self.matrices = selected_matrices
        self.alleles = selected_alleles
    
    @classmethod
    def load(cls, filename, unknown_aa='w', pred_name=None):
        # TODO: very fragile, very propred-specific

        if pred_name is None: pred_name = filename
        
        # Figure out the path to the Rosetta database assuming everything was cloned together
        rosdb = None
        try:
            rosdb = os.path.abspath(os.path.dirname(__file__) + '/../../../database')
            if not os.path.isdir(rosdb): rosdb = None
        except:
            # TODO: print a warning?
            pass
        
        # TODO: look harder for file?
        # Can we point directly to the file?  If so, use it.
        if os.path.exists(filename):
            resolved_fn = filename		
        # If we can automatically find the Rosetta database and the file is in it, use that.
        elif rosdb and os.path.exists(rosdb+'/scoring/score_functions/mhc_epitope/'+filename):
             resolved_fn = rosdb+'/scoring/score_functions/mhc_epitope/'+filename
        # If the "ROSETTA" environment variable is set, try finding the database that way.
        else:
            if not os.getenv('ROSETTA'): raise Exception('Unable to find the database file '+filename+'\n')
            resolved_fn = os.getenv('ROSETTA')+'/main/database/scoring/score_functions/mhc_epitope/'+filename
            if not os.path.exists(resolved_fn): raise Exception('Unable to open file '+filename+', resolved to '+resolved_fn)

        with open(resolved_fn, 'r') as infile:
            def get_line():
                line = '#'
                while len(line)==0 or line[0]=='#':
                    line = infile.readline().strip()
                return line

            # Method
            method = get_line()
            if method != 'propred': raise Exception("ERROR: Unknown epitope predictor " + method);

            # Peptide length
            peptide_length = int(get_line())

            # AA order
            aas = get_line()
            if len(aas) != 20: raise Exception("ERROR: Wrong # AA types for epitope predictor " + aas);

            # Allele-specific info
            nallele = int(get_line())
            matrices = []
            for a in range(nallele):
                # Name
                a_name = get_line()
                # Thresholds
                # Note that the first threshold is for 1%, the second for 2%, etc., so pad with a 0 so that thresholds[i] is the i% threshold
                thresholds = [0] + [float(thresh) for thresh in get_line().split(' ')]
                # PWM
                profile = []
                for p in range(peptide_length):
                    profile.append(dict(zip(aas, [float(w) for w in get_line().split(' ')])))
                # Allele
                matrices.append(AlleleMatrix(a_name, thresholds, profile))
        return cls(pred_name, matrices=matrices, peptide_length=peptide_length, unknown_aa=unknown_aa)
    
    def score_peptide(self, pep):
        try:
            details = [m.score(pep, self.thresh) for m in self.matrices]
            return EpitopeScore(sum(1 for s in details if s is not None), details)
        except KeyError as err: # nancanonical => no epitope
            return EpitopeScore()
    
class Propred (EpitopePredictorMatrix):
    allele_sets = { 
        'all': ['DRB1_0101','DRB1_0301','DRB1_0401','DRB1_0701','DRB1_0801','DRB1_1101','DRB1_1301','DRB1_1501'],
        'southwood98': ['DRB1_0101','DRB1_0301','DRB1_0401','DRB1_0701','DRB1_0801','DRB1_1101','DRB1_1301','DRB1_1501'],
        'test': ['DRB1_0101']
    }
    @classmethod
    def load(cls, unknown_aa='w'):
        return super().load('propred8.txt', unknown_aa=unknown_aa, pred_name='propred8')
