#!/usr/bin/env python

'''Given a database created with setup_screening_project and a screening job file, add ligands for each job based on path information on tagged residues in the database

Author: Sam DeLuca'''

from optparse import OptionParser
import json
from ligand_database import *

def init_options():
    usage = "%prog ligands.db3 input_screening_file.js output_screening_file.js"
    parser=OptionParser(usage)
    return parser
    
    
if __name__ == "__main__":
    
    (options,args) = init_options().parse_args()
    
    database = args[0]
    input_js_file = args[1]
    output_js_file = args[2]
    
    with open(input_js_file) as indata:
        input_job_data = json.load(indata)
    
    
    for params_data in get_params_information(database):
        system = params_data["tag"]
        ligand_file = params_data["params_file"]
        for index,job in enumerate(input_job_data):
            if system == input_job_data[index]["group_name"]:
                ligand_file = ligand_file.replace(".params",".pdb")
                try:
                    input_job_data[index]["ligands"].append(ligand_file)
                except KeyError:
                    input_job_data[index]["ligands"] = [ligand_file]
                
    
    with open(output_js_file,'w') as outdata:
        json.dump(input_job_data,outdata,indent=1)