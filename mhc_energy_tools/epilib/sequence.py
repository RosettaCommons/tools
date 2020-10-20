#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simplistic readers to extract sequence information from various file types.
- Sequences can include '_' characters, treated as noncanonicals
- The headers in fasta-format files can end in "@pos" to start residue numbering there; default 1
- Pretty rudimentary handling of PDB files, padding missing residues with '_' and generally dealing only with the standard twenty 3-letter AA codes

@author: Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu 
"""

import collections, os.path

# The 1-letter codes of the 20 coded AAs
AAs = set('ACDEFGHIKLMNPQRSTVWY')
# A mapping from 3-letter codes to 1-letter codes, with X for unknown
# TODO: include additional 3-1 mappings suitable for epitope prediction (as in CYX below)
AA31 = collections.defaultdict(lambda:'X', 
                               [('ALA','A'),('CYS','C'),('ASP','D'),('GLU','E'),('PHE','F'),('GLY','G'),('HIS','H'),('ILE','I'),('LYS','K'),('LEU','L'),('MET','M'),('ASN','N'),('PRO','P'),('GLN','Q'),('ARG','R'),('SER','S'),('THR','T'),('VAL','V'),('TRP','W'),('TYR','Y')]
                                + [('CYX','C')])

class Sequence (object):
    """A sequence (String) with a name and a chain id, indexed by residue numbers starting at 1 or some other value."""
    def __init__(self, seq, name='[anonymous]', start=1, chain=None, noncanon='silent'):
        self.seq = seq
        self.name = name
        self.start = start
        self.chain = chain
        
        noncanons = [i for i in range(len(seq)) if seq[i] not in AAs]
        if len(noncanons)>0:
            msg = name + ' has noncanonical AA(s): '+','.join('%s%d'%(seq[i],i+start) for i in noncanons)
            if noncanon == 'warn': print(msg)
            elif noncanon == 'error': raise Exception(msg)

    def __getitem__(self, where):
        """Gets the AA at the position / AAs in the slice, accounting for the starting residue number of the sequence."""
        if isinstance(where, slice): 
            return self.seq[(where.start-self.start):(where.stop-self.start):where.step]
        else: 
            return self.seq[where-self.start]
    
    def __len__(self):
        return len(self.seq)
    

def filecore(path):
    """Strips the path to just the core filename, no path, no extension."""
    return os.path.splitext(os.path.split(path)[1])[0]

def extract_name_start(line):
    """From a fasta-format sequence header, gets the name of the sequence. 
    If the name end in '@pos', extract that out as the starting position."""
    name = line[1:]
    at_pos = name.find('@')
    if at_pos >= 0:
        return (name[:at_pos], int(name[at_pos+1:]))
    return (name, 1)
 
def load_fa(filename, noncanon='silent'):
    """Loads modified fasta-format file ('>' optional, '@pos' allowed)"""
    with open(filename, 'rt') as infile:
        rows = infile.readlines()
        if len(rows)==0:
            raise Exception('empty file')
        elif rows[0][0]=='>':
            (name, start) = extract_name_start(rows[0].strip())
            return Sequence(''.join(row.strip() for row in rows[1:]), name=name, start=start, noncanon=noncanon)
        else:
            return Sequence(''.join(row.strip() for row in rows), noncanon=noncanon)
        
def load_pep(filename, noncanon='silent'):
    """Loads a file with a single anonymous peptide per line"""
    seqs = []
    core = filecore(filename)
    with open(filename, 'rt') as infile:
        for (i,seq) in enumerate(infile):
            seqs.append(Sequence(seq.strip(), name=core+'_'+str(i), noncanon=noncanon))
    return seqs

def load_fsa(filename, noncanon='silent'):
    """Loads a modified multiple fasta-format sequence file ('>' required, '@pos' allowed for each)"""
    seqs = []
    name = ''; start = 1; seq = ''
    with open(filename, 'rt') as infile:
        for line in infile:
            line = line.strip()
            if len(line) == 0: continue
            if line[0]=='>':
                if seq != '': seqs.append(Sequence(seq, name=name, start=start, noncanon=noncanon))
                (name,start) = extract_name_start(line)
                seq = ''
            else:
                seq += line
        if seq != '': seqs.append(Sequence(seq, name=name, start=start, noncanon=noncanon))
    return seqs

def load_pdb(filename, noncanon='silent'):
    """Extracts sequence information from CA ATOM records in PDB file.
    Pretty half baked / non-robust.
    Fills missing residues with '_'.
    Separates out different chain ids."""
    seqs = []
    seq = ''
    start = None
    last_chain = None
    last_pos = 0
    core = filecore(filename)
    with open(filename, 'rt') as infile:
        for line in infile:
            if line[:4]=='ATOM' and line[12:16].strip()=='CA':
                chain = line[21]
                pos = int(line[22:26])
                res_type = line[17:20].strip()
                if chain != last_chain:
                    if len(seq)>0:
                        seqs.append(Sequence(seq, name=core+'_'+last_chain, start=start, chain=last_chain, noncanon=noncanon))
                    last_chain = chain
                    seq = ''
                    start = pos
                    last_pos = pos
                elif last_pos is None:
                    last_pos = pos
                elif pos>last_pos+1:
                    seq += '_'*(pos-last_pos-1)
                    if pos==last_pos+2:
                        print("warning: missing amino acid %d in chain %s; treating as epitope-breaking noncanonical" % (last_pos+1, chain))
                    else:
                        print("warning: missing amino acids %d-%d in chain %s; treating as epitope-breaking noncanonical" % (last_pos+1, pos-1, chain))
                    last_pos = pos
                else:
                    last_pos = pos
                seq += AA31[res_type]
    if len(seq)>0: seqs.append(Sequence(seq, name=core+'_'+last_chain, start=start, chain=chain, noncanon=noncanon))
    return seqs
  