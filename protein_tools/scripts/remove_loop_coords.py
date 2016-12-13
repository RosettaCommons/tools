#!/usr/bin/env python

try:
    import rosettautil
except ImportError:
    # if this script is in the Rosetta/tools/protein_tools/scripts/ directory
    # rosettautil is in the ../ directory. Add that to the path. and re-import
    import sys, os
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    import rosettautil

try:
    import Bio.PDB
except ImportError:
    import sys
    sys.stderr.write("\nERROR: This script requires that Biopython (http://biopython.org) is installed.\n\n")
    sys.exit()

import array
from optparse import OptionParser

from rosettautil.protein import util
from rosettautil.rosetta import loops
from rosettautil.util import fileutil

usage = "%prog [options] loopfile.txt input.pdb output.pdb"
parser=OptionParser(usage)
(options,args)=parser.parse_args()

loop_manager = loops.RosettaLoopManager()
loop_manager.read(args[0])
input_struct = util.load_pdb(args[1])

zero_triplet = array.array('f',[0.0,0.0,0.0])

for atom in input_struct.get_atoms():
    resnum = atom.get_parent().get_id()[1]
    if loop_manager.is_res_in_loop(resnum):
        atom.set_coord(zero_triplet)
        atom.set_occupancy(-1.0)

pdb_io = Bio.PDB.PDBIO()
pdb_io.set_structure(input_struct)
outfile = fileutil.universal_open(args[2],'w')
pdb_io.save(outfile)
outfile.close()
