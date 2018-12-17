#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Precomputes epitope scores for a sequence and its considered mutations, storing this information in the database for efficient lookup at design time.

@author: Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu 
"""


# TODO: possible near-term extensions
# - AA choices from resfile
# - fancier scoring functions layered on predictions
# - recompute the score column in a db from the other already-computed columns, with a different scoring function, allele subset/weights, etc.
# - output as csv instead of db -- where to put metadata?
# - merge dbs

import itertools, argparse, csv, functools, operator
from epilib.epitope_predictor_matrix import EpitopePredictorMatrix, Propred
from epilib.epitope_database import EpitopeDatabase
from epilib.netmhcII import NetMHCII
from epilib.sequence import load_fa, load_pdb

def wt_choices(wt):
    """A dictionary of position -> [AA] listing allowed amino acids for each positions, initialized here with just the wild-type at each position"""
    return dict((i, [wt[i]]) for i in range(wt.start, wt.start+len(wt)))

def load_aa_choices_csv(wt, filename):
    """Loads a csv-format file of mutational choices for the wt Sequence. Each row has the position (same numbering as wt) and list of alternative AAs."""
    choices = wt_choices(wt)
    with open(filename, 'r') as infile:
        for row in csv.reader(infile):
            pos = int(row[0]); aas = row[1:]
            if wt[pos] in aas: choices[pos] = aas
            else: choices[pos] += aas
    return choices

def load_aa_choices_pssm(wt, filename, thresh=1):
    """Loads a BLAST PSSM-format file and returns as mutational choices those whose negative log probabilities are at least the threshold."""
    choices = wt_choices(wt)
    with open(filename, 'rt') as infile:
        #Loop until we get to the amino acid order line
        linecount = 0 #Linecounter for AA order line.
        for line in infile:
            #aa_order = "" #aa_order should be reset for each iteration, so that we can identify a failure to find the AA line
            linecount = linecount + 1
            if line == "\n": continue #blank line
            if line.lstrip()[0] == "#": continue #if the first character is a #, treat as a comment
            #The command line psiblast PSSMs have the following header line.  Skip it.
            if line == "Last position-specific scoring matrix computed, weighted observed percentages rounded down, information per position, and relative weight of gapless real matches to pseudocounts": continue
            
            #Split the line into a list of whitespace-delimited strings
            line = line.split()
            
            #The first non-blank, non-comment, non-header line should be the AA line.
            #The command line version of psiblast has the 20 AAs repeated for the two tables in that format
            #The online version has the 20 AAs prefixed by P (position) and C (consensus)
            #Guess the format based on the number of elements in line
            recognized = False
            if len(line) == 40: #Command line psiblast format
                aa_order = line[:20]
                recognized = True
            elif len(line) == 22: #Online psiblast format
                aa_order = line[2:22]
                recognized = True
            elif len(line) == 20: #If it is exaclty 20 elements, let's assume that it's the 20 AAs, with a warning
                aa_order = line
            else:
                aa_order = line
                #If the first 20 or last 20 elements match the 20 AAs, go with it.
                if set(aa_order[:20]) == set('ACDEFGHIKLMNPQRSTVWY'):
                    aa_order = aa_order[:20]
                elif set(aa_order[-20:]) == set('ACDEFGHIKLMNPQRSTVWY'):
                    aa_order = aa_order[-20:]
                else:#If not, continue to the next line
                    continue
            
            if not recognized: print("Line "+str(linecount)+" appears to be the AA order line, but the PSSM is not in a recognized format.  You may get unexpected results!\nSee tools/mhc_energy_tools/pssm_examples for supported formats.")
            
            assert len(aa_order) == 20 #The aa_order must be 20, or something's gone wrong.
            #If the aa_order matches the 20 AAs, break out of this loop and start processing the PSSM contents in the next loop.
            if set(aa_order) == set('ACDEFGHIKLMNPQRSTVWY'): break
        
        #If aa_order not the 20 AAs (i.e. we reached the end of the file without finding it), raise an exception.
        #print('aa order:', aa_order)
        if set(aa_order) != set('ACDEFGHIKLMNPQRSTVWY'):
            raise Exception("Could not process the PSSM format!\n\nSee tools/mhc_energy_tools/pssm_examples for supported formats.")
        
        pos = 0
        for line in infile:
            cols = line.split()
            try:
                pos = int(cols[0])
            except:
                # no longer in matrix; should be finished with sequence
                if pos != wt.start + len(wt) - 1:
                    raise Exception("PSSM stopped at position "+str(pos))
                # finished
                break
            if wt[pos] != cols[1]: raise Exception('mismatched PSSM vs. wt: %s, %s' % (cols[1], wt[pos]))
            scores = [int(cols[j+2]) for j in range(20)]
            aas = [aa_order[j] for j in range(20) if scores[j]>=thresh]
            if wt[pos] in aas: choices[pos] = aas
            else: choices[pos] += aas
    return choices

def generate_peptides(wt, aa_choices, peptide_length, db=None):
    """Given the aa_choices, generates all combinations into peptides of the given length.
    If given a db, only yields those not already there."""
    for pos in range(wt.start, wt.start+len(wt)-peptide_length+1):
        #print('@', pos, '<=', functools.reduce(operator.mul, (len(aa_choices[pos+i]) for i in range(peptide_length))), 'peptides')
        for combo in itertools.product(*[aa_choices[pos+i] for i in range(peptide_length)]):
            peptide = ''.join(combo)
            if db is None or not db.has_peptide(peptide): yield peptide

def setup_parser():
    """Creates an ArgumentParser and sets up all its arguments.
    returns: ArgumentParser"""
    parser = argparse.ArgumentParser(description='Establishes a database of precomputed epitope scores, to enable efficient lookup at design time', 
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog="""Additional notes:
- Currently only handles single sequence per file. Thus if using a pdb file with multiple chains, must specify which with --pdb_chain.
- Pretty rudimentary handling of PDB files, padding missing residues with '_' (i.e., no epitopes) and generally dealing only with the standard twenty 3-letter AA codes
- Wild-type is always considered at each position. If no other choices are provided, only wild-type epitopes are added to the database.
- If the database already exists, it will be augmented. However, the excisting database must be for the same predictor and parameters, including alleles
- NetMHCII can be used in a single run as a subprocess, or in three steps -- use this to save out the peptides (don't provide db), run them through NetMHCII (executable or website), and use this to parse the raw output file
- Res file generation needs chain info, which can be specified by --chain for fasta file and for selection from pdb file
    """)
    # where to get the sequence
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--fa', help="name of file with single protein in fasta-ish format ('>' line optional)")
    source.add_argument('--pdb', help='name of file with single protein embedded in pdb format (allowing chain breaks)')
    parser.add_argument('--chain', help='name of chain to associate with fasta file or to extract from pdb file (when it includes more than one)')
    # positions to mutate
    parser.add_argument('--positions', help='positions to consider mutating, format <start1>[-<end1>],<start2>[-<end2>],... where each comma-separated entry indicates either a single position <start> or an inclusive range <start>-<stop> (default: all positions)')
    # and positions not to mutate (punch holes in the position intervals)
    parser.add_argument('--lock', help='positions not to allow mutating, format <start1>[-<end1>],<start2>[-<end2>],... where each comma-separated entry indicates either a single position <start> or an inclusive range <start>-<stop> (default: none)')
    # where to get the AA choices
    choices = parser.add_mutually_exclusive_group()
    choices.add_argument('--aa_csv', help='name of file with AA choices in csv format; each row has a residue number and list of allowed choices at that position')
    choices.add_argument('--pssm', help='name of BLAST-formatted PSSM')
    # AA choices parameters
    parser.add_argument('--pssm_thresh', help='threshold for AA choices from PSSM: take those with -log prob >= this value (default: %(default)i)', type=int, default=1)
    # epitope predictor
    pred = parser.add_mutually_exclusive_group()
    pred.add_argument('--matrix', help='epitope predictor matrix filename')
    pred.add_argument('--netmhcii', action='store_true', help='use netmhcII executable')
    pred.add_argument('--netmhcii_raw', help='load raw output file from netmhcII corresponding to the input sequence')
    pred.add_argument('--propred', action='store_true', help='use propred matrices')
    # epitope predictor alleles
    alleles = parser.add_mutually_exclusive_group()
    alleles.add_argument('--allele_set', help='name of predefined set of alleles', choices=['test', 'greenbaum11', 'paul15', 'all'])
    alleles.add_argument('--alleles', help='comma-separated list of allele names')
    # epitope predictor parameters
    parser.add_argument('--epi_thresh', help='epitope predictor threshold (default: %(default).2f)', type=int, default=5)
    parser.add_argument('--noncanon', help='how to treat letters other than the 20 canonical AAs (default: %(default)s)', choices=['error', 'silent', 'warn'])
    parser.add_argument('--netmhcii_score', help='type of score to compute (default %(default)s)', choices=['rank','absolute'], default='rank')
    # multiprocessing
    parser.add_argument('--nproc', help='number of processors to distribute peptide predictions across (default: %(default)d)', type=int, default=1)
    parser.add_argument('--batch', help='The number of peptides to distribute to each process for each batch of scoring. When that process finishes, another batch will be sent to a new process until all peptides are scored (default: %(default)d)', type=int, default=10000)
    # optional outputs
    parser.add_argument('--estimate_size', action='store_true', help='print out estimates of numbers of peptides')
    parser.add_argument('--peps_out', help='name of file in which to store raw list of peptide sequences covering choices, with one peptide per line')
    parser.add_argument('--res_out', help='name of file in which to store processed position-specific AA choices, in resfile format; requires chain to be specified by --chain')
    parser.add_argument('--res_header', help='resfile command (e.g. NATRO) to be applied in the resfile header, to be applied to all residues not explicitly specified in the resfile (i.e. with multiple allowed identities).  By default, no command will be included, which is the equivalent to ALLAA.')
    # the database
    parser.add_argument('db', nargs='?', help='name of sqlite3 database to store epitope information (create or augment)')

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
            raise Exception('allele_set '+args.allele_set+' not supported')
        pred.filter_alleles(pred.allele_sets[args.allele_set])
    elif args.alleles is not None:
        pred.filter_alleles(args.alleles.split(','))
    
    # wild-type sequence
    if args.fa is not None:
        wt = load_fa(args.fa)
        wt.chain = args.chain
    elif args.pdb is not None:
        chains = load_pdb(args.pdb)
        if len(chains)==1:
            wt = chains[0]
            if args.chain is not None and args.chain != wt.chain:
                raise Exception('no chain chain in the pdb file')
        else:
            if args.chain is None:
                raise Exception('multi-chain pdb file; please use --chain option to indicate which to use')
            desired_seq = [seq for seq in chains if seq.chain == args.chain]
            if len(desired_seq)==0:
                raise Exception('no chain %s in pdb file' % (args.chain,))
            if len(desired_seq)>1:
                raise Exception('parser fail: multiple chains %s in pdb file' % (args.chain,))
            wt = desired_seq[0]

    # AA choices
    if args.aa_csv is not None:
        aa_choices = load_aa_choices_csv(wt, args.aa_csv)
    elif args.pssm is not None:
        aa_choices = load_aa_choices_pssm(wt, args.pssm, args.pssm_thresh)
    else:
        aa_choices = wt_choices(wt)
        
    # positions
    def get_idx(s):
        try:
            p = int(s)
            assert(p >= wt.start and p <= wt.start+len(wt))
            return p
        except:
            raise Exception('bad position '+s)
    if args.positions is not None:
        # parse which positions can be mutated
        idxs = set()
        for p in args.positions.split(','):
            if '-' in p:
                (start,stop) = p.split('-')
                idxs.update(range(get_idx(start), get_idx(stop)+1))
            else:
                idxs.add(get_idx(p))
        # lock down others
        for i in aa_choices:
            if i not in idxs:
                aa_choices[i] = [wt[i]]
    if args.lock is not None:
        # parse which positions cannot be mutated
        for p in args.lock.split(','):
            if '-' in p:
                (start,stop) = p.split('-')
                for i in range(get_idx(start), get_idx(stop)+1): aa_choices[i] = [wt[i]]
            else:
                i = get_idx(p)
                aa_choices[i] = [wt[i]]

    #print(aa_choices)
    if args.res_header is not None and args.res_out is None:
        raise Exception('--res_out must be specified if you want to apply a --res_header.')
    if args.res_out is not None:
        if wt.chain is None: raise Exception('use --chain to specify chain for res file generation')
        with open(args.res_out, 'w') as outfile:
            # Write the default resfile behaviour, if specified.
            if args.res_header is not None:
                outfile.write(args.res_header + '\n')
            # Start the body of the resfile
            outfile.write('start\n')
            # Loop over positions, and add a PIKAA line for all those with more than 1 AA choice
            # The native AA is always allowed, so if aa_choices[1] == 1, it is just the WT identity
            for i in sorted(aa_choices):
                if len(aa_choices[i])>1: # not just NATAA
                    outfile.write('%d %s PIKAA %s\n' % (i, wt.chain, ''.join(aa_choices[i])))

    if args.estimate_size:
        sizes = dict((pos, functools.reduce(operator.mul, (len(aa_choices[pos+i]) for i in range(pred.peptide_length))))
                     for pos in range(wt.start, wt.start+len(wt)-pred.peptide_length+1))
        print(sizes)
        print('total',sum(sizes.values()))
        if args.db is None and args.peps_out is None:
            # nothing more to do
            return
        
    # open (and maybe initialize) the database
    if args.db is None:
        db = None
    else:
        db = EpitopeDatabase.for_writing(args.db, pred.name, pred.alleles, pred.peptide_length)
        
    if args.netmhcii_raw is not None:
        # process the raw binding predictions and store in the db
        if db is None: raise Exception('need an output database into which to store the raw binding predictions')
        db.save_scores(pred.load_file(args.netmhcii_raw))
        return
    
    if args.peps_out:
        peps_file = open(args.peps_out, 'w')

    # set up multiprocessing if using netmhcii
    if args.nproc > 1 and args.netmhcii is not False:
        print('multiprocessing with',args.nproc,'processors')
        import multiprocessing
        pool = multiprocessing.Pool(processes=args.nproc)
    
    # generate peptides; optional save to file and/or score and save to db
    pep_gen = generate_peptides(wt, aa_choices, pred.peptide_length, None if (db is None or db.is_new) else db)
    while True:
        batches = []
        # If using NetMHCII, get a batch of fresh peptides to send to each proc
        # Otherwise, set the number of processors to 1 and batch size to 0 (don't use batches)
        if args.netmhcii is False:
            nproc = 1
            batch_size = 0
        else:
            nproc = args.nproc
            batch_size = args.batch
        for p in range(nproc):
            # If not using batches, convert pep_gen to a list called batch
            if batch_size == 0:
                batch = list(pep_gen)
            # otherwise, batch should be batch_size peptides from pep_gen
            else:
                batch = list(next(pep_gen) for i in range(batch_size))
            if len(batch)==0: break
            if args.peps_out:
                for peptide in batch: peps_file.write(peptide+'\n')
            batches.append(batch)
        if len(batches)==0: break
        if db is None: continue  # no scoring
        # farm out the scoring
        if nproc == 1:
            results = [pred.score_peptides(batch) for batch in batches]
        else:
            results = pool.map_async(pred.score_peptides, batches).get()
        for scores in results:
            db.save_scores(scores)

    if args.peps_out:
        peps_file.close()

if __name__ == '__main__':
    parser = setup_parser()
    args = parser.parse_args()
    try:
        main(args)
    except Exception as e:
        print('ERROR', e)
        print()
        parser.print_usage()
