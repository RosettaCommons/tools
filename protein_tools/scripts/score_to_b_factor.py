#!/usr/bin/env python2.7
'''
score_to_b_factor.py -- Convert the score table in a Rosetta output to B-factors.

From the Meiler Lab script repository.
'''

import sys
from optparse import OptionParser
from Bio.PDB import *
from Bio.PDB import PDBExceptions
import warnings
warnings.simplefilter('ignore',PDBExceptions.PDBConstructionWarning)
from multiprocessing import Pool as poolio

try:
    import rosettautil
except ImportError:
    # if this script is in the tools/protein_tools/scripts/ directory
    # rosettautil is in the ../ directory. Add that to the path. and re-import
    import sys, os
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    import rosettautil
import rosettautil.rosetta
print rosettautil.rosetta.__file__
import rosettautil.rosetta.rosettaScore_beta

usage = "%prog <pdbs>  use --help for more options"
parser = OptionParser(usage)
parser.add_option("--term",dest="term",default="total",help="defaults to total the scoring term to set b factors to")
parser.add_option("--score_list",dest="score_list",help="instead of passing rosetta score term, pass list of residues with scores you want mapped",default="")
parser.add_option("--keep-table",dest="table",help="Preserve the rosetta score table at the bottom of the pdb",action="store_true",default=False)
parser.add_option("--prefix",dest="prefix",help="the prefix on the output file, defaults to B_factor_",default="b_factor_")
parser.add_option("--log",dest="log",help="If your ranges between highest score and lowest score are very wide, you can plot them as a log base where -log <base> is the log base you want to use",type=int)
parser.add_option("--multi",dest="multi_proc",help="Multiprocess the job.",action="store_true", default=False)
(options,args)=parser.parse_args()

print "WARNING: it is absolutely imperative that the residue numbers in your pdb file start with 1 and proceed with no gaps!"

if len(args) < 1:
    parser.error("specify input pdb(s)")

io=PDBIO()

def strip_path(pdb):
    return pdb.split('/')[-1]

def score_to_bfactor(pdb):
    scores = rosettaScore.ScoreTable(pdb)
    #print "got scores"
    PDBParse = PDBParser(PERMISSIVE=1)
    struct = PDBParse.get_structure(pdb[0:3],pdb)
    print "got structures"
    atoms = struct.get_atoms()
    print "filling atoms"
    map={}
    if(options.score_list!=""):
        file=open(options.score_list,'r')
        for line in file:
            fields=line.split()
            residue=int(fields[0])
            access=float(fields[1])
            map[residue]=access
            file.close()
    for atom in atoms:
        #print "parsing atoms"
        residue_id = atom.get_parent().get_id()[1]
        chain_id = atom.get_parent().get_parent().get_id()
        if(options.score_list!=""):
            try:
                score=map[residue_id]
                print "setting b_factor"
                atom.set_bfactor(score)
            except KeyError:
                atom.set_bfactor(0.0)
        else:
            score = scores.get_score(chain=chain_id,pdbres=residue_id,term=options.term)[0]
            if options.log and score != 0:
                import math
                if score < 0:
                    atom.set_bfactor(-math.log(math.fabs(score),options.log))
                    print "set a negative b fact"
                else:
                    atom.set_bfactor(math.log(score,options.log))
                    print "set a negative b fact"
            else:
                atom.set_bfactor(score)
    io.set_structure(struct)
    io.save(options.prefix+strip_path(pdb))
    if(options.table):
        raw_table = rosettaScore.get_table(pdb)
        outfile = open(options.prefix+strip_path(pdb),'a')
        outfile.writelines(raw_table)
        outfile.close()

if __name__ == "__main__":
    if options.multi_proc:
        try:
             poolio().map(score_to_bfactor, args)
        except KeyboardInterrupt:
            p.terminate()
            quit()
    else:
        for i in args:
            score_to_bfactor(i)
