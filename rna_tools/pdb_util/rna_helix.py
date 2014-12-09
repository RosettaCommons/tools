#! /usr/bin/env python
import argparse
import tempfile
import subprocess
import shutil
from renumber_pdb_in_place import renumber_pdb
from parse_options import get_resnum_chain
from rosetta_exe import rosetta_exe

parser = argparse.ArgumentParser(description='Run rna_helix to build a helix')
parser.add_argument('-seq', required=True, help='Sequence of the helix, in the format: ggaa uucc', nargs=2 )
parser.add_argument('-resnum', nargs='+',help='Renumber the residues with input format: 13-16 30-33')
parser.add_argument('-o', default='helix.pdb', help='Filename of output pdb')
parser.add_argument('-weights', default='', help='Weights file defining score function')
parser.add_argument('-finish_weights', default='', help='Weights file defining a finisher score function')
parser.add_argument('-silent', default='', help='silent file output')
parser.add_argument('-dump', action='store_true', default=False, help='dump intermediate pdbs')
parser.add_argument('-put_intra_into_total', action='store_true',default=False, help='calculate intra-res terms and include in totals')
args = parser.parse_args()


#Build the helix
temp = tempfile.NamedTemporaryFile(delete=False)
cmdline  = rosetta_exe('rna_helix')
cmdline += (' -rna::corrected_geo  '+
            '-score:rna_torsion_potential RNA11_based_new ' +
            '-chemical::enlarge_H_lj ')
cmdline += '-o %s ' % temp.name
cmdline += '-seq '
for i in args.seq:    cmdline += i + ' '
if len( args.weights ) > 0:        cmdline += '-score:weights %s ' % args.weights
if len( args.finish_weights ) > 0: cmdline += '-finish_weights %s ' % args.finish_weights
if len( args.silent ) > 0: cmdline += '-out:file:silent %s ' % args.silent
if args.dump: cmdline += '-dump '
if args.put_intra_into_total: cmdline += '-put_intra_into_total '

print 'Rosetta cmdline:', cmdline
subprocess.check_call(cmdline.split())

output_file = args.o
# slightly weird

#Renumber the pdb
if args.resnum is not None:
    resnums = []
    chains  = []
    for i in args.resnum:  get_resnum_chain(i, resnums, chains )
    renumber_pdb([temp.name], resnums, chains)

temp.close()
shutil.move(temp.name, output_file)

