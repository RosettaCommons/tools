#!/usr/bin/python3

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


###############################################################################
### helper functions
###############################################################################
def safe_submit(command, verbose = False):
    if verbose:
    	print(command if isinstance(command, str) else ' '.join(command))
    out, err = sp.Popen(
        command,
        shell=isinstance(command, str),
        stdout=sp.PIPE,
        stderr=sp.PIPE
    ).communicate()
    if err:
        print('\n'.join([out, err]))
        return False
    return out


###############################################################################
### main functions
###############################################################################
def rna_thread(options):
    options = init_options(options)
    rna_thread = ["rna_thread"]
    rna_thread.append('-s %s' % options.in_file_s)
    rna_thread.append('-o %s' % options.out_file_o)
    if options.fasta:
        rna_thread.append('-fasta %s' % options.in_file_fasta)
    if options.seq:
        rna_thread.append('-seq %s' % options.seq)
    if options.seq_offset:
        rna_thread.append('-seq_offset %s' % options.seq_offset)
    if options.residue_type_set:
        rna_thread.append('-in:file:residue_type_set %s' % options.residue_type_set)
    if options.extra_res:
        rna_thread.append('-in:file:extra_res %s' % options.extra_res)
    if options.extra_res_fa:
        rna_thread.append('-in:file:extra_res_fa %s' % options.extra_res_fa)    
    if options.options:
        rna_thread += ['-' + o for o in options.options.split('-') if len(o)]
    out = safe_submit(rna_thread, options.verbose)
    if options.verbose:
        print(out)
    if not os.path.exists(options.out_file_o):
        return None
    return options.out_file_o 


###############################################################################
### init functions
###############################################################################
def init_options(options):
    if isinstance(options, list):
        parser = init_options_parser()
        options = parser.parse_args(args = options)
    if options.in_file_s is None:
        print("[WARNING] Must supply a template pdb (-s) for threading!")
        sys.exit(1)
    if not len([_f for _f in [options.fasta, options.seq] if _f]):
        print("[WARNING] Must specify either -fasta, -seq or -mutate")
        sys.exit(1)
    return options

def init_options_parser():
    parser = argparse.ArgumentParser(
        description='Wrapper for rna_thread.cc',
        add_help=False
    )
    parser.add_argument(
        "in_file_s",
        help="Template pose '-s'",
    )
    parser.add_argument(
        "-o","-out:file:o",
        dest="out_file_o",
        help="Outfile",
        default="threaded.pdb"
    )
    parser.add_argument(
        "-fasta","-in:file:fasta",
        dest="fasta",
        help="FASTA file",
        default=None
    )
    parser.add_argument( 
        "-seq",
        help="target sequence (can include dashes). Must specify either this or -fasta. Length of this sequence must exactly equal length of template pose specified by -s",
        default=None
    )
    parser.add_argument(
        "-seq_offset",
        help="Integer to add to all residue numbers in output PDB",
        default=None
    )
    parser.add_argument(
        "-sequence_mask_file",
        help="Output name for sequence mask file (not in use at the moment)",
        default=None
    )
    parser.add_argument(
        "--residue_type_set", '-in:file:residue_type_set',
        dest="residue_type_set",
        help="residue type set, used in rosetta",
        default=None
    )
    parser.add_argument(
        "--extra_res", '-in:file:extra_res',
        dest="extra_res",
        help="location of .params file for a residue type, used in rosetta",
        default=None
    )
    parser.add_argument(
        "--extra_res_fa", '-in:file:extra_res_fa',
        dest="extra_res_fa",
        help="location of .params file for a (full-atom) residue type, used in rosetta",
        default=None
    )
    parser.add_argument(
        "-options",
        help="Extra Rosetta options",
        default=None
    )
    parser.add_argument(
        "-v","--verbose",
        help="Increase verbosity",
        action='store_true'
    )
    return parser

###############################################################################
### main
###############################################################################
if __name__=='__main__':
    
    print(rna_thread(sys.argv[1:]))
