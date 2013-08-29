#!/usr/bin/env python

import sys
import json
import glob
import re
from optparse import OptionParser

def init_options():
    usage = "%prog params_dir/ structure_dir/ output.js"
    parser=OptionParser(usage)
    return parser

if __name__ == "__main__":
    
    options,args = init_options().parse_args()
    params_dir = args[0]
    structure_dir = args[1]
    output_js_path = args[2]

    pdb_map = {}
    for path in glob.glob(params_dir+"/*/*.params"):
        base_name = path.split("/")[-1].split(".")[0]
        base_path = path.split("/")
        base_path.pop()
        pdb_path = "/".join(base_path)+"/"+base_name+".pdb"
        with open(path) as params_file:
            for line in params_file:
                line = line.split()
                if len(line) >= 3:
                    if line[1] == "system_name":
                        try:
                            pdb_map[line[2]].append(pdb_path)
                        except KeyError:
                            pdb_map[line[2]] = [pdb_path]
    jobs = []
    total_jobs = 0
    for system_name in pdb_map:
        new_job = {"group_name" : system_name, "ligands" : pdb_map[system_name]}
        protein_list = []
        for path in glob.glob(structure_dir+"/*"+system_name+"_*.pdb*"):
            protein_list.append(path)
    
        new_job["proteins"] = protein_list
        total_jobs += len(pdb_map[system_name])*len(protein_list)
    
        jobs.append(new_job)
    
    with open(output_js_path,'w') as output_file:
        json.dump(jobs,output_file,indent=1)
        
    print total_jobs,"Jobs"