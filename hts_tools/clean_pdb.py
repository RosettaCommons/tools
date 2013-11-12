#!/usr/bin/env python2.7

from Bio.PDB import * 
from optparse import OptionParser
import sys
import json 
import os

class CleanSelect(Select):
    def __init__(self,chains,ligands):
        self.chains = chains
        self.ligands = ligands
    def accept_residue(self,residue):
        resname = residue.get_resname()
        resname = resname.strip()
        if resname in Polypeptide.aa3 or resname in self.ligands:
            return 1
        else:
            return 0
                
    def accept_chain(self,chain):
        if chain.get_id() in self.chains:
            return 1
        else:
            return 0
        

def center_residue(residue):
    coordinate_sum = [0.0,0.0,0.0]
    coordinate_count = 0.0
    for atom in residue:
        coordinate_sum += atom.get_coord()
        coordinate_count += 1.0
    
    return coordinate_sum/coordinate_count

usage = "%prog --chains=A,B --ligands_to_keep=CL,ABC --ligand_for_center=LIG input.pdb output.pdb center_file.js"
parser = OptionParser(usage)
parser.add_option("--chains",dest="chains",default="")
parser.add_option("--ligands_to_keep",dest="ligands",default="")
parser.add_option("--ligand_for_center",dest="center_ligands",default="")
(options,args) = parser.parse_args()

input_path = args[0]
output_path = args[1]
center_file = args[2]

chains = options.chains.split(",")
ligands = options.ligands.split(",")

if os.path.exists(center_file):
    with open(center_file) as infile:
        center_map = json.load(infile)
else:
    center_map = {}

base_name = input_path.split("/")[-1].split(".")[0]

parser = PDBParser()
structure = parser.get_structure("input",input_path)

for residue in structure.get_residues():
    chain = residue.get_full_id()[2]
    residue_name = residue.get_resname()
    if chain in chains and residue_name == options.center_ligands:
        center_map[base_name] = center_residue(residue).tolist()

with open(center_file,'w') as outfile:
    json.dump(center_map,outfile,indent=1)

io=PDBIO()
io.set_structure(structure)
io.save(output_path,CleanSelect(chains,ligands))


