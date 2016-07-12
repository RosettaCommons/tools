#!/usr/bin/env python2.7
'''USAGE: %prog [options] infile.pdb outfile.cst

Generate a Rosetta constraint file to restrain atom pair distances.

By default, will generate a constraint file with HARMONIC constraints for all CA pairs which fall within the specified cutoff.

If generating constraints for docking purposes, the --partners option can specify the partners such that constraints will not be generated across the interface.

Additionally, individual residues can be left as unconstrained (through --release_file), or can be selected to have a different constraint strength (with --relax_file).

*IMPORTANT* This script assumes the input PDB has been renumbered such that the PDB numbering matches the Pose numbering. *IMPORTANT*
'''

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

from optparse import OptionParser
import warnings
from rosettautil.protein import util
from rosettautil.util import fileutil
from math import log
import subprocess
import sys

################## options for the numerable datapoints required for this analysis ##############################################

parser=OptionParser(__doc__)
parser.add_option("--distance_cutoff","--distance","--cutoff",dest="dist_cutoff",
    help="cutoff distance between residue Cbeta, default: %default",default=8.0,type=float)
parser.add_option("--standard_deviation","--stdev",dest="stdev",
    help="the standard deviation for the HARMONIC constraint, default: %default",default=0.5,type=float)
parser.add_option("--release_docking_partners", "--partners", dest="partners",
    help="Do not constrain across the interface. Chains of each partner are specfied in the Rosetta style (AB_CD)",default="")
parser.add_option("--release_file",dest="release_file",
    help="file with column of integers for residues which should not be restrained", default="")
parser.add_option("--relax_file",dest="relax_file",
    help="file with column of integers for residues which should have relaxed restraints", default="")
parser.add_option("--relax_value",dest="relax_value",
    help="standard deviation for the HARMONIC constraint on residues in the relax file. default: %default", default=0.5,type=float)
parser.add_option("--constrain_atom",dest="constrain_atom",
    help="which atom to constrain, default: %default",default="CA")

(options,args) = parser.parse_args()
warnings.simplefilter('ignore',Bio.PDB.PDBExceptions.PDBConstructionWarning)

if len(args) != 2:
    parser.error("You must specify both the reference PDB structure and the name of the output constraint file.")

structure = util.load_pdb(args[0])                #read in the pdb to calculate residue neigbors

# Only apply constraints to the atoms with the specified name.
structure_atoms = [atom for atom in structure.get_atoms() if atom.get_id() == options.constrain_atom ]

if len(structure_atoms) < 2:
    parser.error("There are fewer than two '"+options.constrain_atom+"' atoms in the structure - can't generate constraints.")

partner1 = []
partner2 = []
if not options.partners == "":
    partner_groups = options.partners.split('_')
    if len(partner_groups) != 2:
        parser.error("The --partners option must have one (and only one) underscore character.")
    for chain_id in partner_groups[0]:
        partner1.append(chain_id)
    for chain_id in partner_groups[1]:
        partner2.append(chain_id)

relax_list = []
if not options.relax_file == "":
    for line in open(options.relax_file, 'r'):
        relax_list.append(line.strip())
release_list = []
if not options.release_file == "":
    for line in open(options.release_file, 'r'):
        release_list.append(line.strip())

#enumerate provides both an index and a list
#this will let us iterate through only half the array, instead of
#processing every constraint twice
outfile = open(args[1],'w')
for index, atom in enumerate(structure_atoms):                    # iterate across atoms first time
    res_num = atom.get_full_id()[3][1]
    for atom2 in structure_atoms[index:]:                    # iterate across atoms that have not already been evaluated relative to first iteration
        atom_chain = atom.get_parent().get_parent().get_id()    # get parent chain of atom in first iteration loop
        atom2_chain = atom2.get_parent().get_parent().get_id()    # get parent chain of atom2 in second iteration loop
        res_num2 = atom2.get_full_id()[3][1]
        if res_num == res_num2:                            # Do not restrain an atom to itself. That's just silly.
            continue
        if str(res_num) in release_list or str(res_num2) in release_list :
            continue
        elif ( (atom_chain in partner1) and (atom2_chain in partner2) ) or ( (atom_chain in partner2) and (atom2_chain in partner1) ):    # check to see if atoms are across the interface from each other. If true then do not apply a restraint.
            continue
        else:
            dist = atom2 - atom            # calculate distance between atoms
            if dist <= options.dist_cutoff and (str(res_num) in relax_list or str(res_num2) in relax_list):
                outfile.write(str("AtomPair "+atom.get_id()+" "+str(res_num)+" "+atom2.get_id()+" "+str(res_num2)+" HARMONIC "+str(round(dist,3))+" "+str(options.relax_value)+"\n"))
            elif dist <= options.dist_cutoff:        # only apply restraints to atoms within cutoff distance
                outfile.write(str("AtomPair "+atom.get_id()+" "+str(res_num)+" "+atom2.get_id()+" "+str(res_num2)+" HARMONIC "+str(round(dist,3))+" "+str(options.stdev)+"\n"))        # write the atom_pair_constraint line
