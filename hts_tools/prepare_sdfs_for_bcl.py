#!/usr/bin/env python

'''Given a database created with setup_screening_project and a json file containing rosetta descriptor data, produce sdf files with Rosetta descriptors formatted for input into the BCL

Author: Sam DeLuca'''

from optparse import OptionParser
import json
from sdf_parser import SdfFile
from ligand_database import *
import os

def init_options():
    usage = "%prog ligands.db3 descriptor_data.js output.sdf"
    parser=OptionParser(usage)
    parser.add_option("--clobber",dest="clobber",default=False, action="store_true")
    return parser
    
if __name__ == "__main__":
    parser = init_options()
    options, args = parser.parse_args()
    
    if len(args) != 3:
        parser.exit("You must specify an hts pipeline database, a Rosetta descriptor JSON file and a path for an output pdf")
        
    database_path = args[0]
    descriptor_path = args[1]
    output_path = args[2]
    
    if os.path.exists(output_path) and not options.clobber:
        parser.exit(output_path+" already exists, use --clobber to overwrite")
    elif os.path.exists(output_path) and options.clobber:
        os.remove(output_path)
    
    with open(descriptor_path) as descriptor_file:
        descriptor_data = json.load(descriptor_file)
        
    for record in descriptor_data:
        ligand_name = record["name3"]
        sdf_file = get_sdf_path_from_ligand_name(database_path,ligand_name)
        mol_data = SdfFile.parse_sdf(sdf_file).next() #just get the first conformer
        
        for key in record:
            value = [record[key]] #interface to molfile needs a list
            mol_data.add_data_entry(key,value)
        
        
        with open(output_path,'a') as output_file:
            for line in mol_data.data_to_string_list():
                output_file.write(str(line)+"\n")
            output_file.write("$$$$\n")