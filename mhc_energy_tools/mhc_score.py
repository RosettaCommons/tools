#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Performs basic epitope prediction for a sequence or set of sequences, to support pre- and post-processing of mhc_epitope_energy optimized designs.

@author: Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu 
"""

import argparse, csv, sys, os.path
from epilib.epitope_predictor_matrix import EpitopePredictorMatrix, Propred
from epilib.epitope_csv import EpitopeCSV
from epilib.epitope_database import EpitopeDatabase
from epilib.netmhcII import NetMHCII
from epilib.sequence import Sequence, load_fa, load_pep, load_fsa, load_pdb
from epilib.ext_citations import CitationTracker

def setup_parser():
    """Creates an ArgumentParser and sets up all its arguments.
    returns: ArgumentParser"""
    parser = argparse.ArgumentParser(description="""Predicts epitopes for a sequence or set of sequences. 
Performs basic epitope prediction for a sequence or set of sequences, to support pre- and post-processing of mhc_epitope energy optimized designs.
Provide sequence(s) and specify predictor and its parameters; the script generates either the total score (useful for testing) or a peptide-by-peptide report or plot.""", 
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog="""Additional notes:
- Only uses one source for sequences, checking in order: command-line, specified input file, stdin
- Sequences can include '_' characters, treated as noncanonicals, which in the current implementation forces the containing peptides to be non-epitopes
- The "fasta-ish" header can end in "@pos" to start residue numbering there; default 1
- Pretty rudimentary handling of PDB files, padding missing residues with '_' (i.e., no epitopes) and generally dealing only with the standard twenty 3-letter AA codes
- If no predictor specified, defaults to propred
- Allele sets include 'all' and 'test', specific to predictors, as well as two published lists that are supported by NetMHCII:
  + greenbaum11: www.ncbi.nlm.nih.gov/pubmed/21305276
  + paul15: www.ncbi.nlm.nih.gov/pubmed/25862607
- Looks for predictor matrix file in current directory as well as in $ROSETTA/main/database/scoring/score_terms/mhc_epitope.
- Output filenames can include a '$' character that is substituted with the current sequence's name, for cases where the input includes multiple sequences""")
    # where to get the sequence(s), unless on command line
    source = parser.add_mutually_exclusive_group()
    source.add_argument('--fa', help="name of file with single sequence in fasta-ish format, '>' line optional")
    source.add_argument('--fsa', help="name of file with one or more sequences in fasta-ish format, each including '>' line")
    source.add_argument('--pdb', help='name of file with one or more sequences embedded in pdb format')
    source.add_argument('--pep', help="name of file with one sequence per line")
    parser.add_argument('--chain', help='if only care about one chain and pdb file has more than one, the id of the desired one')
    # epitope predictor
    pred = parser.add_mutually_exclusive_group()
    pred.add_argument('--csv', help='name of csv database from which to load epitope predictions')
    pred.add_argument('--db', help='name of database from which to load epitope predictions')
    pred.add_argument('--matrix', help='generic epitope predictor matrix filename')
    pred.add_argument('--netmhcii', action='store_true', help='use netmhcII executable')
    pred.add_argument('--propred', action='store_true', help='use propred matrices')
    # epitope predictor alleles
    alleles = parser.add_mutually_exclusive_group()
    alleles.add_argument('--allele_set', help='name of predefined set of alleles', choices=['test', 'greenbaum11', 'paul15', 'southwood98', 'all'])
    alleles.add_argument('--alleles', help='comma-separated list of allele names')
    # epitope predictor parameters
    parser.add_argument('--epi_thresh', help='epitope predictor threshold (default: %(default).2f)', type=int, default=5)
    parser.add_argument('--noncanon', help='how to treat letters other than the 20 canonical AAs (default: %(default)s)', choices=['error', 'silent', 'warn'])
    parser.add_argument('--netmhcii_score', help='type of score to compute (default %(default)s)', choices=['rank','absolute'], default='rank')
    parser.add_argument('--db_unseen', help='how to handle unseen epitope (default %(default)s)', choices=['warn','error','score'], default='warn')
    parser.add_argument('--db_unseen_score', help='what score to use for unseen epitope (default %(default)i)', type=int, default=100)
    # output
    parser.add_argument('--report', help='level of detail for output report: just _total_ score / epitope _hits_ / scores for _all_ peptides (default: %(default)s)', choices=['total','hits','full'], default='total')
    parser.add_argument('--report_file', help='filename in which to save report in csv format (default: print to stdout); see note about wildcard')
    parser.add_argument('--plot_hits_file', help='filename in which to save plot of epitope hits (format specified by extension); see note about wildcard')
    # sequence(s) on command-line
    parser.add_argument('seq', help='raw sequence as string', nargs='*')
    
    return parser

def handle_seq(seq, pred, args):
    """Helper to process a single Sequence."""
    epimap = pred.score_protein(seq.seq)
    print('*', seq.name, epimap.total_score())
    # Report
    if args.report=='full' or args.report=='hits':
        # One row per peptide (full or just hits), starting with a header
        rows = [['pos','peptide','score'] + epimap.alleles]
        total = 0
        for (i, peptide) in enumerate(epimap.peptides):
            score = epimap.p2s[peptide]
            total += score.value
            if args.report=='full':
                rows.append([i+seq.start, peptide, score.value] + [str(h) if h is not None else '-' for h in score.details])
            elif score.value>0:
                rows.append([i+seq.start, peptide, score.value] + [str(h) if (h is not None and h<=pred.thresh) else '-' for h in score.details])
        # Save or print
        if args.report_file is not None:
            fn = args.report_file.replace('$', seq.name)
            with open(fn, 'a') as outfile:
                outfile.write(seq.name + '\n')
                csv.writer(outfile, lineterminator=os.linesep).writerows(rows)
                outfile.write('\n')
            print('  saved to '+fn)
        else:
            csv.writer(sys.stdout).writerows(rows)
    # Plot
    if args.plot_hits_file is not None:
        # TODO: allow customization of size, color, labels, ....
        try:
            import matplotlib.pyplot as plt
        except:
            print('matplotlib is required for plotting, and could not be imported.  Make sure it is correctly installed.\n')
            exit()
        plt.clf()
        plt.bar(range(seq.start, seq.start+len(epimap.peptides)), [epimap.peptide_score(p).value for p in epimap.peptides])
        fn = args.plot_hits_file.replace('$', seq.name)
        # Check if the fn already exists
        if os.path.exists(fn):
            print('WARNING: Attempting to save an epitope hit plot for sequence ' + fn + ' to a location with an existing file (' + fn + ').  The plot will not be generated.')
            print('WARNING: You may have forgotten to include a "$" in your filename if you are scoring multiple sequences.')
        # If not, generate the plot
        else:
            plt.savefig(fn)
            print('  plotted to '+fn)
        
def main(args):
    """Sets up an epitope prediction run based on the parsed args.
    args: argparse.Namespace"""
    
    # Citation manager object to track what features are being used
    cite_manager = CitationTracker()
    
    # epitope predictor
    if args.matrix is not None:
        pred = EpitopePredictorMatrix.load(args.matrix)
        #We don't have any generic matrix strategy.  If added, add the citation.
    elif args.netmhcii:
        cite_manager.add_citation(cite_manager.netmhcii)
        pred = NetMHCII(score_type=args.netmhcii_score[0])
    elif args.csv is not None:
        #Generic, so we can't add a citation
        pred = EpitopeCSV.for_reading(args.csv, handle_unseen=args.db_unseen[0], unseen_score=args.db_unseen_score)
    elif args.db is not None:
        #Generic, so we can't add a citation
        pred = EpitopeDatabase.for_reading(args.db, handle_unseen=args.db_unseen[0], unseen_score=args.db_unseen_score)
    else:
        cite_manager.add_citation(cite_manager.propred)
        pred = Propred.load()
    pred.thresh = args.epi_thresh
    
    # alleles
    if args.allele_set is not None:
        if args.allele_set not in pred.allele_sets: 
            raise Exception('allele_set '+args.allele_set+' not supported')
        pred.filter_alleles(pred.allele_sets[args.allele_set])
    elif args.alleles is not None:
        pred.filter_alleles(args.alleles.split(','))

    # Get sequences from whatever source is specified
    # Predict epitopes for each
    if len(args.seq)>0:
        for (i, seq) in enumerate(args.seq): handle_seq(Sequence(seq, name='seq'+str(i)), pred, args)           
    elif args.fa is not None:
        handle_seq(load_fa(args.fa), pred, args)
    elif args.fsa is not None:
        for seq in load_fsa(args.fsa): handle_seq(seq, pred, args)
    elif args.pdb is not None:
        got_one = False
        for seq in load_pdb(args.pdb): 
            if args.chain is None or seq.chain==args.chain:
                got_one = True
                handle_seq(seq, pred, args)
        if not got_one: raise Exception('no chain %s in pdb file' % (args.chain,))
    elif args.pep is not None:
        for seq in load_pep(args.pep): handle_seq(seq, pred, args)
    else:
        print('enter sequences:')
        i = 0
        for line in sys.stdin:
            handle_seq(Sequence(line.strip(), name='seq'+str(i)), pred, args)
            i += 1
    
    cite_manager.output_citations()

if __name__ == '__main__':
    parser = setup_parser()
    args = parser.parse_args()
    try:
        main(args)
    except Exception as e:
        print('ERROR', e)
        print()
        parser.print_usage()  
