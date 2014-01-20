#!/usr/bin/env python2.7
import json
from optparse import OptionParser
import glob
import math 
import sys
import os
from Bio import PDB
import numpy

def get_emptiest_bin(bin_map):
    smallest_bin = 0
    smallest_bin_size = 1000000000
    for bin_id in bin_map:
        bin_size = 0
        for record in bin_map[bin_id]:
            bin_size += record[0]
        if smallest_bin_size > bin_size:
            smallest_bin = bin_id
            smallest_bin_size = bin_size
    return smallest_bin

def find_ligand_center(pdb_path):
    parser = PDB.PDBParser()
    structure = parser.get_structure("x",pdb_path)
    atom_sum = numpy.zeros(3)
    atom_count = 1.0
    for atom in structure.get_atoms():
        atom_sum += atom.get_coord()
        atom_count += 1.0
    return atom_sum/atom_count

def parse_params_dir(params_dir,group_name,group_list=None):
    pdb_map = {}
    params_map = {}
    found_param = False
    for path in glob.glob(params_dir+"/*/*.params"):
        found_param = True
        base_name = path.split("/")[-1].split(".")[0]
        base_path = path.split("/")
        base_path.pop()
        pdb_path = "/".join(base_path)+"/"+base_name+".pdb"
        if group_name == None:
            with open(path) as params_file:
                for line in params_file:
                    line = line.split()
                    if len(line) >= 3:
                        if line[1] == "system_name":
                            if (group_list != None and line[2] in group_list) or group_list == None:
                                try:
                                    pdb_map[line[2]].append(pdb_path)
                                except KeyError:
                                    pdb_map[line[2]] = [pdb_path]
                            
        else:
            try:
                pdb_map[group_name].append(pdb_path)
            except KeyError:
                pdb_map[group_name] = [pdb_path]
                            

        params_map[base_name] = path
    if found_param == False:
        sys.exit("ERROR: %s contains no params files" % params_dir)
    return (pdb_map,params_map)
    
    
def chunks(l, n):
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def make_inactive_crossdock_jobs(pdb_map,structure_dir,cutoff):
    jobs = []
    total_jobs = 0
    for system_name in pdb_map:
        protein_list = []
        short_system_name = system_name.split("_")[0]
        for path in glob.glob(structure_dir+"/*"+short_system_name+"*.pdb*"):
            protein_list.append(path)
        
        ligands_to_dock = []
        for cross_name in pdb_map:
            if cross_name == system_name:
                ligand_center = find_ligand_center(pdb_map[system_name][0])
            ligands_to_dock += pdb_map[cross_name]
        if len(protein_list) == 0 or len(ligands_to_dock) == 0:
            continue
        ligands_per_job = int(math.ceil(cutoff/len(protein_list)))
        for ligand_chunk in chunks(ligands_to_dock,ligands_per_job):
            new_job = {"group_name" : system_name}
            new_job["proteins"] = protein_list
            new_job["ligands"] = ligand_chunk
            new_job["startfrom"] = ligand_center.tolist()
            jobs.append(new_job)
    return jobs
        
def make_jobs(pdb_map,structure_dir,cutoff,group_name,native_dir =None):

    jobs = []
    total_jobs = 0
    for system_name in pdb_map:
        protein_list = []
        if group_name == None:
            #short_system_name = system_name.split("_")[0]
            for path in glob.glob(structure_dir+"/*"+system_name+"*.pdb*"):
                protein_list.append(path)
        else:
            short_system_name = group_name
            for path in glob.glob(structure_dir+"/*.pdb*"):
                protein_list.append(path)
        if len(protein_list) == 0 or len(pdb_map[system_name]) == 0:
            continue
        ligands_per_job = int(math.ceil(cutoff/len(protein_list)))
        for ligand_chunk in chunks(pdb_map[system_name],ligands_per_job):
            new_job = {
                "group_name" : system_name,
                "proteins" : protein_list,
                "ligands" : ligand_chunk
            }
            if native_dir != None:
                if len(ligand_chunk) == 0:
                    sys.exit("ERROR: the option --create_native_commands doesn't work if a job has more than 1 ligand associated with it")
                ligand_path = ligand_chunk[0]
                for path in glob.glob(native_dir+"/"+system_name+"*.pdb*"):
                    native_pdb = path
                    break
                new_job["native"] = "%s %s" % (native_pdb,ligand_path)
            jobs.append(new_job)
    return jobs

def get_params_for_job_list(params_map,job_list):
    params_list = []
    for job in job_list:
        base_names = [ x.split("/")[-1].split(".")[0] for x in job["ligands"] ]
        for name in base_names:
            params_list.append(params_map[name])
    return params_list
    

def init_options():
    usage = "%prog --n_chunks param_dir structure_dir output_prefix"
    parser=OptionParser(usage)
    parser.add_option("--n_chunks",dest="n_chunks",help="number of total files to split work into. (Default: 4)",default=4)
    parser.add_option("--max_per_job",dest="max_per_job",help="maximum number of structures per job. (Default: 10000)",default=10000)
    parser.add_option("--group_name",dest="group_name",help="Manually specify a group name to be used for all proteins and all ligands")
    parser.add_option("--group_name_list",dest="group_name_list",help="specify a list of acceptable group names seperated by commas, jobs will only be constructed from those groups",default=None)
    parser.add_option("--inactive_cross_dock",dest="cross_dock",help="create jobs that dock all actives as inactives",default=False, action="store_true")
    parser.add_option("--create_native_commands",dest="native_dir",help="path to a directory of pdb files containing native proteins. (Optional)",default=None)
    return parser
    

if __name__ == "__main__":
    
    (options,args) = init_options().parse_args()
    param_dir = args[0]
    structure_dir = args[1]
    prefix = args[2]
    bin_count = int(options.n_chunks)
    cutoff = int(options.max_per_job)
    
    if options.group_name_list != None:
        if os.path.exists(options.group_name_list):
            group_name_list = []
            with open(options.group_name_list) as infile:
                for line in infile:
                    line = line.rstrip()
                    group_name_list.append(line)
        else:
            group_name_list = options.group_name_list.split(",")
    else:
        group_name_list = []
    
    
    if options.group_name_list == None:
        pdb_map,params_map = parse_params_dir(param_dir,options.group_name)
    else:
        pdb_map,params_map = parse_params_dir(param_dir,options.group_name,group_list=group_name_list)
    if options.cross_dock:
        jobs = make_inactive_crossdock_jobs(pdb_map,structure_dir,cutoff)
    else:
        jobs = make_jobs(pdb_map,structure_dir,cutoff,options.group_name,native_dir=options.native_dir)
    

    job_list = []
    all_jobs = 0
    for job in jobs:
        n_proteins = len(job["proteins"])
        n_ligands = len(job["ligands"])
        total_size = n_proteins*n_ligands
        job_list.append( (total_size,job) )
        all_jobs += total_size

    job_list = sorted(job_list,key=lambda x: x[0],reverse=True)

    bin_data = {}
    for bin_id in xrange(1,bin_count+1):
        bin_data[bin_id] = []

    for job_size,job_data in job_list:
        #print job_size, job_data
        bin_to_use = get_emptiest_bin(bin_data)
        bin_data[bin_to_use].append((job_size,job_data))

    print "Job files have the following sizes:"
    for bin_id in bin_data:
        bin_size = 0
        for record in bin_data[bin_id]:
            bin_size += record[0]
        print bin_id,bin_size
            
    for bin_id in bin_data:
        new_jobs = [record[1] for record in bin_data[bin_id]]
        params = get_params_for_job_list(params_map,new_jobs)
        complete_file_object = {
            "params" : params,
            "jobs" : new_jobs
        }
        with open("%s_%02d.js" % (prefix,bin_id),'w') as outfile:
            json.dump(complete_file_object,outfile,indent=1)
