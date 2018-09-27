#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Precomputes epitope scores for a sequence and its considered mutations, storing this information in the database for efficient lookup at design time.
Provide a sequence and specification of possible mutations; the script generates and scores epitopes and creates/augments the database.
Currently supports prediction via
- matrix methods, with code particularly targeted to Propred matrices (as provided with the distribution),
- the NetMHCII executable, wrapping the invocation and processing the results

@author: Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu 
"""


# TODO: possible near-term extensions
# - multiple chains in a single file
# - AA choices from resfile
# - other PSSM formats (though I think this is the one supported in Rosetta itself)
# - fancier scoring functions layered on predictions
# - recompute the score column in a db from the other already-computed columns, with a different scoring function, allele subset/weights, etc.

import itertools, argparse, csv
from epilib.epitope_predictor_matrix import EpitopePredictorMatrix, Propred
from epilib.epitope_database import EpitopeDatabase
from epilib.netmhcII import NetMHCII
from epilib.sequence import load_fa, load_pdb

def generate_peptides(aa_choices, peptide_length):
    """Given the aa_choices [ [AA] ], generates all combinations into peptides of the given length."""
    peptides = set()
    for i in range(len(aa_choices)-peptide_length+1):
        aa_combos = list(itertools.product(*[aa_choices[i+p] for p in range(peptide_length)]))
        peptides.update(''.join(combo) for combo in aa_combos)
    return peptides

def save_peptides(peptides, filename):
    """Writes the peptides to the file."""
    with open(filename, 'wt') as of:
        for pep in peptides:
            of.write(pep)
            of.write('\n')

def load_aa_choices_csv(wt, filename):
    """Load a csv-format file of mutational choices for the wt Sequence. Each row has the position (same numbering as wt) and list of alternative AAs."""
    choices = [[aa] for aa in wt] # default unless specifed
    with open(filename, 'r') as infile:
        for row in csv.reader(infile):
            pos = int(row[0]); aas = row[1:]
            if wt[pos] in aas: choices[pos] = aas
            else: choices[pos] += aas
    return choices

def load_aa_choices_pssm(wt, filename, thresh=0):
    """Loads a BLAST PSSM-format file and returns as mutational choices those whose negative log probabilities are at least the threshold."""
    choices = [[aa] for aa in wt] # default unless specifed
    with open(filename, 'rt') as infile:
        infile.readline() # blank
        infile.readline() # header
        line = infile.readline() # AA order
        aa_order = line.split()[:20]
        for line in infile:
            cols = line.split()
            pos = int(cols[0])-1
            if wt[pos] != cols[1]: raise Exception('mismatched PSSM vs. wt: %s, %s' % (cols[1], wt[pos]))
            scores = [int(cols[j+2]) for j in range(20)]
            aas = [aa_order[j] for j in range(20) if scores[j]>=thresh]
            if wt[pos] in aas: choices[pos] = aas
            else: choices[pos] += aas
    return choices

def setup_parser():
    """Creates an ArgumentParser and sets up all its arguments.
    returns: ArgumentParser"""
    parser = argparse.ArgumentParser(description='Establishes a database of precomputed epitope scores, to enable efficient lookup at design time', 
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog="""Additional notes:
- Currently only handles single sequence per file.
- Wild-type is always considered at each position. If no other choices are provided, only wild-type epitopes are added to the database.
- If the database already exists, it must be for the same predictor and parameters, including alleles
- NetMHCII can be used in a single run as a subprocess, or in three steps -- use this to save out the peptides, run them through (executable or website), and use this to parse the raw output file
    """)
    # where to get the sequence
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--fa', help="name of file with single protein in fasta-ish format ('>' line optional)")
    source.add_argument('--pdb', help='name of file with single protein embedded in pdb format (allowing chain breaks)')
    # where to get the AA choices
    choices = parser.add_mutually_exclusive_group()
    choices.add_argument('--aa_csv', help='name of file with AA choices in csv format; each row has a residue number and list of allowed choices at that position')
    choices.add_argument('--pssm', help='name of BLAST-formatted PSSM')
    # TODO choices.add_argument('--resfile', help='name of file with AA choices in resfile format')
    # AA choices parameters
    parser.add_argument('--pssm_thresh', help='threshold for AA choices from PSSM: take those with -log prob >= this value (default: %(default).2f)', type=float, default=0)
    parser.add_argument('--peps_out', help='name of file in which to store raw list of peptide sequences covering choices, with one peptide per line')
    # epitope predictor
    pred = parser.add_mutually_exclusive_group()
    pred.add_argument('--matrix', help='epitope predictor matrix filename')
    pred.add_argument('--netmhcii', action='store_true', help='use netmhcII executable')
    pred.add_argument('--netmhcii_raw', help='load raw output file from netmhcII corresponding to the input sequence')
    pred.add_argument('--propred', action='store_true', help='use propred matrices')
    # epitope predictor alleles
    alleles = parser.add_mutually_exclusive_group()
    alleles.add_argument('--allele_set', help='name of predefined set of alleles', choices=['test', 'greenbaum11', 'paul15'])
    alleles.add_argument('--alleles', help='comma-separated list of allele names')
    # epitope predictor parameters
    parser.add_argument('--epi_thresh', help='epitope predictor threshold (default: %(default).2f)', type=int, default=5)
    parser.add_argument('--noncanon', help='how to treat letters other than the 20 canonical AAs (default: %(default)s)', choices=['error', 'silent', 'warn'])
    parser.add_argument('--netmhcii_score', help='type of score to compute (default %(default)s)', choices=['rank','absolute'], default='rank')
    parser.add_argument('db', help='name of sqlite3 database to store epitope information (create or augment)')

    return parser

def main(args):
    """Sets up an epitope prediction run based on the parsed args.
    args: ArgumentParser"""
    
    # epitope predictor
    if args.matrix is not None:
        pred = EpitopePredictorMatrix.load(args.matrix)
    elif args.netmhcii or args.netmhcii_raw is not None:
        pred = NetMHCII(score_type=args.netmhcii_score[0])
    elif args.netmhcii_raw is not None:
        pred = NetMHCII(nm_bin=False, score_type=args.netmhcii_score[0])
    else:
        pred = Propred.load()
    pred.thresh = args.epi_thresh

    # alleles
    if args.allele_set is not None:
        if args.allele_set not in pred.allele_sets: 
            raise Exception('alleleset '+args.allele_set+' not supported')
        pred.filter_alleles(pred.allele_sets[args.allele_set])
    elif args.alleles is not None:
        pred.filter_alleles(args.alleles.split(','))
    
    # wild-type sequence
    if args.fa is not None:
        wt = load_fa(args.fa)
    elif args.pdb is not None:
       wt = load_pdb(args.pdb)

    # AA choices
    if args.aa_csv is not None:
        aa_choices = load_aa_choices_csv(wt, args.aa_csv)
    elif args.pssm is not None:
        aa_choices = load_aa_choices_pssm(wt, args.pssm, args.pssm_thresh)
    else:
        aa_choices = [[aa] for aa in wt]
        
    # open (and maybe initialize) the database
    db = EpitopeDatabase.for_writing(args.db, pred.name, pred.alleles, pred.peptide_length)
    
    if args.netmhcii_raw is not None:
        # process the raw binding data and store it in the db
        db.save_scores(pred.load_file(args.netmhcii_raw))
    else:
        # generate and score peptides, and store in the db
        peptides = generate_peptides(aa_choices, pred.peptide_length)
        if args.peps_out: save_peptides(peptides, args.peps_out)
        db.save_scores(pred.score_peptides(peptides))
    
if __name__ == '__main__':
    parser = setup_parser()
    args = parser.parse_args()
    try:
        main(args)
    except Exception as e:
        print('ERROR', e)
        print()
        parser.print_usage()
