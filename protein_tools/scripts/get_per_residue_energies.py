#!/usr/bin/env python2.7

import sys, os
if sys.version_info < (2, 7):
    raise Exception("You must use python2.7 to run this")

import warnings
import argparse
import collections
from multiprocessing import Pool
warnings.simplefilter('ignore', PDBExceptions.PDBConstructionWarning)

try:
    import rosettautil
except ImportError:
    # if this script is in the Rosetta/tools/protein_tools/scripts/ directory
    # rosettautil is in the ../ directory. Add that to the path. and re-import
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import rosettautil

from rosettautil.protein import util, pdbStat

try:
    from rosettautil.protein import amino_acids
except ImportError:
    #Okay, do we have one in the same directory as this script?
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    import amino_acids

try:
    from Bio.PDB import PDBExceptions
except ImportError:
    import sys
    sys.stderr.write("\nERROR: This script requires that Biopython (http://biopython.org) is installed.\n\n")
    sys.exit()


# usage
usage = "%prog [options] *pdbs"
parser = argparse.ArgumentParser(prog="get_per_residue_energies.py",
                                 description="Parse per residue energies into a tab seperated file that is easy to parse")
parser.add_argument("pdbs", metavar="*pdb",nargs=argparse.REMAINDER, help="The pdb files to parse")
parser.add_argument("--multi", "-m", dest="multi", help="run using multiple processors", action='store_true', default=False)
parser.add_argument("--out", "-o", dest="output", help="the name of the  output file", default="per_residue_energies.csv")
parser.add_argument("--tag", "-t", dest="tag",help="A comma seperated keyvalue pair you want to put in the file")
args = parser.parse_args()

def main(pdb):
    pdb_loaded = util.load_pdb(pdb)
    residues = [i for i in pdb_loaded.get_residues()]
    scores = collections.OrderedDict()
    open_pdb = open(pdb)
    score_titles = ""
    breaker = False
    counter = 0
    chain_numbering = 1
    for line in open_pdb:
        try:
		line.split()[0]
	except:
		continue
	if line.split()[0] == 'label':
            score_titles = line.split()[1:]
        if line.split()[0] == 'pose':
            breaker = True
            continue
        if breaker:
            if line.split()[0] == "#END_POSE_ENERGIES_TABLE":
                breaker = False
                continue
            if ':' in line.split()[0]:
                res_id = line.split()[0].split(':')[0]
            else:
                res_id = line.split()[0].split('_')[0]

            try:
                one_letter = amino_acids.longer_names[res_id]
            except KeyError:
                print "Warning - Can't find one letter code for {0}".format(res_id)
                one_letter = res_id

            scores[line.split()[0]] = {'id':line.split()[0],
                                       'res_id':res_id,
                                       'res_id_one_letter':one_letter,
                                       'chain':'',
                                       'chain number': '',
                                       'pose':line.split()[0].split("_")[-1],
                                       'scores':line.split()[1:]}
            if "CtermProteinFull" in line.split()[0]:
                counter += 1
                chain_numbering = 1
    for res_id,identifier in zip(residues,scores):
        scores[identifier]['chain'] = res_id.get_parent().get_id()
	if res_id.get_id()[2].strip():
		scores[identifier]['chain number'] = str(res_id.get_id()[1]) + res_id.get_id()[2]
	else:
		scores[identifier]['chain number'] = res_id.get_id()[1]
    return [pdb.split('/')[-1],scores,score_titles]

def format_and_return(scores,outfile):
    titles = scores[0][2]
    if args.tag:
        tag_title = args.tag.split(",")[0]
        tag = args.tag.split(",")[1]
        with open(outfile,'w') as f:
            header = ['pdb','full_id','res_id','res_one','chain','pose_number','res_number',tag_title] + titles
            f.write(",".join(header))
            f.write("\n")
            for entry in scores:
                scores_dict = entry[1].values()
                for residue in scores_dict:
                    f.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(
                        entry[0],residue['id'],residue['res_id'],
                        residue['res_id_one_letter'],residue['chain'],
                        residue['pose'],residue['chain number'],tag,",".join(residue['scores'])))
    else:
        with open(outfile,'w') as f:
            header = ['pdb','full_id','res_id','res_one','chain','pose_number','res_number'] + titles
            f.write(",".join(header))
            f.write("\n")
            for entry in scores:
                scores_dict = entry[1].values()
                for residue in scores_dict:
                    f.write("{0},{1},{2},{3},{4},{5},{6},{7}\n".format(
                        entry[0],residue['id'],residue['res_id'],
                        residue['res_id_one_letter'],residue['chain'],
                        residue['pose'],residue['chain number'],",".join(residue['scores'])))


if __name__ == '__main__':
    if args.multi:
        try:
            p = Pool()
            scores = p.map(main, args.pdbs)
        except KeyboardInterrupt:
            p.terminate()
        format_and_return(scores,args.output)
    else:
        scores = map(main,args.pdbs)
        format_and_return(scores,args.output)
