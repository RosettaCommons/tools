#!/usr/bin/env python

from optparse import OptionParser
import sys
import json 
import os
import glob

usage = "%prog startfrom_map.js input_proteins/ output_startfrom.js"
parser = OptionParser(usage)
#parser.add_option("--chains",dest="chains",default="")
#parser.add_option("--ligands_to_keep",dest="ligands",default="")
#parser.add_option("--ligand_for_center",dest="center_ligands",default="")
(options,args) = parser.parse_args()

startfrom_map_path = args[0]
protein_dir = args[1]
output_startfrom_path = args[2]

with open(startfrom_map_path) as startfrom_map_file:
    startfrom_map = json.load(startfrom_map_file)

startfrom_list = []
for protein_path in glob.glob(protein_dir+"/*"):
    base_name = protein_path.split("/")[-1].split(".")[0]
    tag = base_name.split("_")[0]
    
    coords = startfrom_map[tag]
    
    #i was stupid and made the startfrom map output a
    #different format than the rosetta startfrom file
    
    startfrom_object = {
        "x" : coords[0],
        "y" : coords[1],
        "z" : coords[2],
        "file_name" : base_name
    }
    
    startfrom_list.append(startfrom_object)
    
with open(output_startfrom_path, 'w') as output_startfrom_file:
    json.dump(startfrom_list,output_startfrom_file,indent=1)