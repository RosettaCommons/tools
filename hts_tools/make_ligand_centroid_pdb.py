#!/usr/bin/env python 

'''
Given a mol2 file make a pdb file taht just has a VRT residue at the geometric center of the residue.  Eventually add mol support too.  This is a preprocessing step for making startfrom files easily.

Author:Sam DeLuca
'''

from optparse import OptionParser


def get_xyz_from_mol2(mol2_path):
    '''def really primitive mol2 parser that just gets atom triples mol2 files'''
    with open(mol2_path) as mol2_file:
        in_atom_block = False
        for line in mol2_file:
            fields = line.split()
            if len(fields) < 1:
                continue
            
            if in_atom_block and len(fields) >= 5:
                yield (float(fields[2]),float(fields[3]),float(fields[4]))
            
            if fields[0] == "@<TRIPOS>ATOM":
                in_atom_block = True
            elif fields[0] == "@<TRIPOS>BOND":
                break
                
                
def init_options():
    usage = "%prog ligand.mol2 outfile.pdb"
    parser=OptionParser(usage)
    return parser
    
if __name__ == "__main__":
    parser = init_options()
    (options,args) = parser.parse_args()
    
    if len(args) != 2:
        parser.error("You must specify both an input file and an output file")
    count = 0.0
    x_total = 0.0
    y_total = 0.0
    z_total = 0.0
    for x,y,z in get_xyz_from_mol2(args[0]):
        x_total += x
        y_total += y
        z_total += z
        count += 1.0
        
    virt_string = "HETATM%5i %-4.4s %3.3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2.2s  \n" % (1, "ORIG", "VRT", "X", 1, x_total/count, y_total/count, z_total/count, 1.0, 20.0, "")

    with open(args[1],'w') as outfile:
        outfile.write(virt_string+"\n")
        
        