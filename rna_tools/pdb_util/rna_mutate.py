#!/usr/bin/python

###############################################################################
### imports
###############################################################################
import sys
import subprocess as sp
import os
import os.path
import shutil
import glob
import argparse
import rna_thread
import get_sequence
import parse_tag


###############################################################################
### helper functions
###############################################################################
def mutate_sequence(residues, mutations, pdb):
    sequences, chains, resnums = get_sequence.get_sequences(pdb) 
    for residue, mutation in zip(residues,mutations):    
        for seqidx, sequence in enumerate(sequences):
            sequence = bytearray(sequence)
            for residx, res in enumerate(sequence):
                if residue[0] != resnums[seqidx][residx]:
                    continue
                if residue[1] != '' and residue[1] != chains[seqidx][residx]:
                    continue
                sequence[residx] = mutation
            sequences[seqidx] = str(sequence)
    return ''.join(sequences)

###############################################################################
### main functions
###############################################################################
def rna_mutate(options):
    options = init_options(options)
    residues = zip(*parse_tag.parse_tag(options.residues))
    mutations = list(options.mutations.lower())
    options.seq = mutate_sequence(residues, mutations, options.in_file_s)
    rna_thread.rna_thread(options)
    return options.out_file_o 

###############################################################################
### init functions
###############################################################################
def init_options(options):
    if isinstance(options, list):
        parser = init_options_parser()
        options = parser.parse_args(args = options)
    return options

def init_options_parser():
    parser = argparse.ArgumentParser(
        parents=[rna_thread.init_options_parser()]
    )
    parser.add_argument(
        "-r","--residues",
        help="Residues to mutate.",
    )
    parser.add_argument(
        "-m","--mutations",
        help="New residue types.",
    )
    return parser
###############################################################################
### main
###############################################################################
if __name__=='__main__':
    
    print rna_mutate(sys.argv[1:])
