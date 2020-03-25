#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written for NetMHCII version 2.3
Can integrate with local binary or just load files saved from the website.
To integrate with the binary, either pass in the path to it or set the environment variable NMHOME should be set to the same value as in the installed netMHC-2.3 shell script

@author: Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu 
"""

import collections, os, subprocess, tempfile
from epilib.epitope_predictor import EpitopePredictor, EpitopeScore, EpitopeMap

class NetMHCII (EpitopePredictor):
    version = '2.3'

    # names as NetMHC wants them -- note HLA for DP/DQ but not for DR
    allele_sets = {
        'test':
            ['DRB1_0101'],
            
        'all':
            ['DRB1_0101','DRB1_0103','DRB1_0301','DRB1_0401','DRB1_0402','DRB1_0403','DRB1_0404','DRB1_0405','DRB1_0701','DRB1_0801','DRB1_0802','DRB1_0901','DRB1_1001','DRB1_1101','DRB1_1201','DRB1_1301','DRB1_1302','DRB1_1501','DRB1_1602',
             'DRB3_0101','DRB3_0202','DRB3_0301','DRB4_0101','DRB4_0103','DRB5_0101',
             'H-2-IAb','H-2-IAd','H-2-IAk','H-2-IAs','H-2-IAu','H-2-IEd','H-2-IEk',
             'HLA-DPA10103-DPB10201','HLA-DPA10103-DPB10301','HLA-DPA10103-DPB10401','HLA-DPA10103-DPB10402','HLA-DPA10103-DPB10601','HLA-DPA10201-DPB10101','HLA-DPA10201-DPB10501','HLA-DPA10201-DPB11401','HLA-DPA10301-DPB10402',
             'HLA-DQA10101-DQB10501','HLA-DQA10102-DQB10501','HLA-DQA10102-DQB10502','HLA-DQA10102-DQB10602','HLA-DQA10103-DQB10603','HLA-DQA10104-DQB10503','HLA-DQA10201-DQB10202','HLA-DQA10201-DQB10301','HLA-DQA10201-DQB10303','HLA-DQA10201-DQB10402','HLA-DQA10301-DQB10301','HLA-DQA10301-DQB10302','HLA-DQA10303-DQB10402','HLA-DQA10401-DQB10402','HLA-DQA10501-DQB10201','HLA-DQA10501-DQB10301','HLA-DQA10501-DQB10302','HLA-DQA10501-DQB10303','HLA-DQA10501-DQB10402','HLA-DQA10601-DQB10402']
,
        # Greenbaum J, Sidney J, Chung J, Brander C, Peters B, Sette A.
        # Functional classification of class II human leukocyte antigen (HLA) molecules reveals seven different supertypes and a surprising degree of repertoire sharing across supertypes
        # Immunogenetics. 2011 Jun;63(6):325-35
        # https://www.ncbi.nlm.nih.gov/pubmed/21305276
        'greenbaum11':
            ['DRB1_0101', 'DRB1_0301', 'DRB1_0401', 'DRB1_0405', 'DRB1_0701', 'DRB1_0802', 'DRB1_0901', 'DRB1_1101', 'DRB1_1201', 'DRB1_1302', 'DRB1_1501',
             'DRB3_0101', 'DRB3_0202', 'DRB4_0101', 'DRB5_0101',
             'HLA-DPA10103-DPB10201', 'HLA-DPA10103-DPB10401', 'HLA-DPA10201-DPB10101', 'HLA-DPA10201-DPB10501', 'HLA-DPA10201-DPB11401', 'HLA-DPA10301-DPB10402',
             'HLA-DQA10101-DQB10501', 'HLA-DQA10102-DQB10602', 'HLA-DQA10301-DQB10302', 'HLA-DQA10401-DQB10402', 'HLA-DQA10501-DQB10201', 'HLA-DQA10501-DQB10301'],

        # Paul S, Lindestam Arlehamn CS, Scriba TJ, Dillon MB, Oseroff C, Hinz D, McKinney DM, Carrasco Pro S, Sidney J, Peters B, Sette A.
        # Development and validation of a broad scheme for prediction of HLA class II restricted T cell epitopes
        # J Immunol Methods. 2015 Jul;422:28-34
        # https://www.ncbi.nlm.nih.gov/pubmed/25862607
        'paul15':
            ['DRB1_0301', 'DRB1_0701', 'DRB1_1501', 'DRB3_0101', 'DRB3_0202', 'DRB4_0101', 'DRB5_0101'],
        
        # Southwood S, Sidney J, Kondo A, del Guercio MF, Appella E, Hoffman S, Kubo RT, Chesnut RW, Grey HM, Sette A.
        # Several common HLA-DR types share largely overlapping peptide binding repertoires
        # J Immunol. 1998 Apr;160:3363-73
        # https://www.ncbi.nlm.nih.gov/pubmed/9531296
        'southwood98':
            ['DRB1_0101','DRB1_0301','DRB1_0401','DRB1_0701','DRB1_0801','DRB1_1101','DRB1_1301','DRB1_1501']
        }

    nm2std = {} # cache for std_name
    @staticmethod
    def std_name(nm_name):
        """Converts NetMHC name into something suitable for a column name in the database.
        - drops "HLA-"
        - converts "-" to "_"
        Caches the conversion."""
        if nm_name in NetMHCII.nm2std: return NetMHCII.nm2std[nm_name]
        name = nm_name.replace('HLA-','').replace('-','_')
        NetMHCII.nm2std[nm_name] = name
        return name
    
    def __init__(self, nm_alleles=None, thresh=5, score_type='r', nm_bin=None):
        """nm_alleles uses NetMHCII naming, as in the predefined sets above
        score_type 'r' means relative (percentile) and 'a' means absolute ('IC50'); note that a different thresh would be appropriate for 'a'
        nm_bin is the executable; else looks for it in $NMHOME"""
        if nm_alleles is None: nm_alleles = NetMHCII.allele_sets['paul15']
        super().__init__('netmhcii+'+NetMHCII.version, alleles=[NetMHCII.std_name(a) for a in nm_alleles], peptide_length=15)
        self.nm_alleles = nm_alleles
        self.thresh = thresh
        self.score_type = score_type
        if nm_bin is not None:
            nm_bin = nm_bin
        else:
            path = os.getenv('NMHOME')
            if path is None: raise Exception('please set the NMHOME environemnt variable to the path where the netMHC-2.3 executable lives')
            self.nm_bin = path+'/netMHCII-2.3'
 
    def filter_alleles(self, selected_alleles):
        missing = set(selected_alleles) - set(NetMHCII.allele_sets['all'])
        if len(missing)>0:
            raise Exception('no such allele(s) '+missing)
        self.nm_alleles = selected_alleles
        self.alleles = [NetMHCII.std_name(a) for a in selected_alleles]

    def process_raw(self, rows, is_peptides):
        """Processes the raw output from NetMHCII.
        is_peptides indicates whether it was run in peptide mode (a single peptide per line) or protein mode (a whole sequence in a fasta file)"""
        got_version = False
        peptides = []
        scores = collections.defaultdict(dict)
        # Set which columns contain the allele, position, peptide sequence, and scores
        # These column numbers may change if NetMHCII output changes in a future version
        allele_col = 0 #Allele is always the first column
        pos_col = 1 #Position is always the second column
        pep_col = 2 #Peptide is always the third column
        #Score column changes depending on whether scoring peptides or a protein sequence, and whether using relative or absolute scores
        #If using absolute scoring, the score column is 5 (log-transformed: 1-log50k(aff))
        if self.score_type == 'a': 
            score_col = 5
        #If using relative scoring with peptides, the score column is 8
        elif is_peptides:
            score_col = 8
        #If using relative scoring with a full protein sequence, the score column is 7
        else:
            score_col = 7
        for row in rows:
            if len(row)==0 or row[0]=='#': continue
            cols = row.split()
            if len(cols)==0: continue
            if 'version' in cols:
                # try to make sure it's the same version
                if cols[-1].startswith(self.version):
                    got_version = True
                else:
                    print('*** untested version ',cols[-1],'; use at your own risk!')
            elif cols[allele_col] in self.nm_alleles:
                # a score row for an allele -- extract and store the info if good enough
                score = float(cols[score_col])
                #if score <= self.thresh:
                allele = NetMHCII.std_name(cols[allele_col])
                peptide = cols[pep_col]
                idx = int(cols[pos_col]) - 1
                #The idx basically gives the position number of the current line in the sequence or in the list of peptides.
                #If the number of peptides is greater than the index, we should no longer be on the first allele
                #That means that the peptide has already been included in the peptides list
                if idx < len(peptides):
                    # Since we are on a non-first allele, the peptide should already be present in the list at that index.
                    # If not, raise an exception.
                    if peptide != peptides[idx]:
                        raise Exception('peptide mismatch for %d: %s vs. %s' % (idx, peptides[idx], peptide))
                #If the idx number and the length of peptides is the same, we are on the first allele and the peptide hasn't been seen.
                #Add it to peptides.
                elif idx == len(peptides):
                    peptides.append(peptide)
                #If the idx number is greater than the length of peptides, that means we are missing peptides for some reason.
                #Raise an exception.
                else:
                    raise Exception('missing peptides from %d to %d' % (len(peptides), idx))
                scores[peptide][allele] = score
        if not got_version: print('*** unidentified version, use at your own risk!')
        episcores = []
        for peptide in peptides:
            # NetMHC converts non-20 AAs to 'X' and gladly scores them as contributing to epitope score
            # Override to non-epitope
            # TODO: in peptide mode, don't waste effort scoring
            if 'X' in peptide: 
                episcores.append(EpitopeScore())
                continue
            details = [] #Store the details of all peptide/allele combos
            meet_thresh = 0 #Store the number of peptide/allele combos that are less than self.thresh
            for allele in self.alleles:
                details.append(scores[peptide][allele])
                if scores[peptide][allele] <= self.thresh:
                    meet_thresh+=1
            episcores.append(EpitopeScore(meet_thresh, details))
        return EpitopeMap(self.peptide_length, self.alleles, peptides, episcores)

    def load_file(self, filename):
        """Saved-out NetMHCII file."""
        with open(filename, 'r') as infile:
            return self.process_raw(infile.readlines())

    def run(self, filename, is_peptides=True):  
        """Runs the executable as a subprocess."""
        cmd = [self.nm_bin] + (['-p'] if is_peptides else []) + ['-a', ','.join(self.nm_alleles), '-f', filename]
        print('invoking NetMHCII: ' + ' '.join(cmd))
        sp = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        (nm_stdout, nm_stderr) = sp.communicate()
        if nm_stderr != '': print('*** netmhc stderr:', nm_stderr)
        #print('*** nm stdout:', nm_stdout)
        #TODO: save nm output if desired???
        return self.process_raw(nm_stdout.split('\n'), is_peptides)
        
    def score_peptide(self, pep):
        return self.score_peptides([pep])
    
    def score_peptides(self, peptides, filename=None):
        """Stores the peptides in a temporary file unless a filename is given; invokes the executable; collects the scores."""
        is_tmp = filename is None
        if is_tmp:
            with tempfile.NamedTemporaryFile(delete=False, mode='wt') as fp:
                filename = fp.name
                print('saving peptides in tempfile', fp.name)
                for peptide in peptides:
                    fp.write(peptide)
                    fp.write('\n')
        epimap = self.run(filename, is_peptides=True)
        if is_tmp: os.unlink(filename)
        return epimap
    
    def score_protein(self, seq, filename=None):
        """Stores the protein in a temporary file unless a filename is given; invokes the executable; collects the scores."""
        is_tmp = filename is None
        if is_tmp:
            with tempfile.NamedTemporaryFile(delete=False, mode='wt') as fp:
                filename = fp.name
                print('saving sequence in tempfile', fp.name)
                fp.write('>tmp\n')
                fp.write(seq)
                fp.write('\n')
        epimap = self.run(filename, is_peptides=False)
        if is_tmp: os.unlink(filename)
        return epimap
