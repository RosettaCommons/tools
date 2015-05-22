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
def mutate_sequence(residue, mutation, pdb):
    sequences, chains, resnums = get_sequence.get_sequences(pdb) 
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
    residue = tuple([x[0] for x in parse_tag.parse_tag(options.residue)])
    mutation = options.mutation.lower()
    options.seq = mutate_sequence(residue, mutation, options.in_file_s)
    if "z" in options.seq:
        options.residue_type_set = "rna_ligands"
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
        "residue",
        help="Residue to mutate.",
    )
    parser.add_argument(
        "mutation",
        help="New residue type.",
    )
    return parser
###############################################################################
### main
###############################################################################
if __name__=='__main__':
    
    print rna_mutate(sys.argv[1:])
