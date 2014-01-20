#!/usr/bin/env python

'''Given a database created with setup_screening_project, make params files for every ligand. Params files will be stored in a filename hashed directory structure.

Author: Sam DeLuca'''

from math import ceil, log
import itertools 
import os
import hashlib
import subprocess
from multiprocessing import Pool
from optparse import OptionParser
from ligand_database import * 
import shutil
path_to_params = ""

def get_name_from_params(path,database):
    '''Given the path to a params file, return the IO_STRING value.  if it doesn't exist, return None'''
    path=path.strip()
    #if path is commented, ignore it
    if len(path) == 0:
        return 0
    if path[0] == "#":
        return 0
    fullpath = database+path
    #print fullpath
    path = path.split("/")
    #is path a params file?
    type = path[-1].split(".").pop()
    if type == "params":
        if os.path.exists(fullpath):
            paramsfile = open(fullpath,'r')
            for line in paramsfile:
                line = line.split()
                if len(line) >0:
                    if line[0] == "NAME":
                        return line[1]
        else:
            return None
    else:
        return None

def get_disallowed_ligands(database):
    '''Return a set of 3 letter names which are already assigned to ligands in the database'''
    residue_type_set_path = database+"/chemical/residue_type_sets/"
    residue_set_list = os.listdir(residue_type_set_path)
    disallowed_ligands = set()
    for residue_set in residue_set_list:
        if residue_set[0] != ".":
            current_residue_set_path = residue_type_set_path+residue_set+"/"
            residue_types = open(current_residue_set_path+"residue_types.txt",'r')
            for line in residue_types:
                name = get_name_from_params(line, current_residue_set_path)
                #print name
                if name != None:
                    disallowed_ligands.add(name)
    return disallowed_ligands


def required_ligand_name_length(ligand_count):
    '''Get the necessary length of param file name given the number of ligands being processed'''
    character_space = 36 #[A-Z0-9]
    return max(2,int(ceil(log(ligand_count)/log(character_space))))
 
def setup_dir_for_filehash(output_dir):
    '''Setup all the directories used in filehashing ahead of time so i dont have to worry about my threads blocking properly'''
    #i bet theres a better way to do this next line, whatevers
    character_set = ['a','b','c','d','e','f','0','1','2','3','4','5','6','7','8','9']
    subdirs = ["".join(y) for y in itertools.product(character_set,repeat=2)]
    for subdir in subdirs:
        outpath = output_dir+"/"+subdir
        if not os.path.exists(outpath):
            os.mkdir(outpath)
            
def make_param_file(input_sdf_path,ligand_name,output_dir):
    sha1 = hashlib.sha1()
    sha1.update(ligand_name)
    digest = sha1.hexdigest()
    sub_dir = digest[0:2]
    full_dir = output_dir+"/"+sub_dir
    param_file_command = path_to_params + " -n %(name)s --long-names --clobber --conformers-in-one-file --mm-as-virt %(path)s"
    
    output_param_path = "%s.params" % ligand_name
    output_pdb_path = "%s.pdb" % ligand_name
    output_conformer_path = "%s_conformers.pdb" % ligand_name
    
    subprocess.call(param_file_command % {"name" : ligand_name, "path" : input_sdf_path},shell=True)
    
    try:
        shutil.move(os.getcwd()+"/"+output_param_path,full_dir + "/" + output_param_path)
    except IOError:
        print "something went wrong with",ligand_name
        return None
    shutil.move(os.getcwd()+"/"+output_pdb_path,full_dir + "/" + output_pdb_path)
    try:
        shutil.move(os.getcwd()+"/"+output_conformer_path,full_dir + "/" + output_conformer_path)
    except IOError:
        print "No conformers for",ligand_name
    return full_dir + "/" + output_param_path
    
    
def append_activity_tag_to_params(data):
    if data["output_path"] != None:
        with open(data["output_path"],'a') as paramfile:
            paramfile.write("STRING_PROPERTY %s %s\n" % ("system_name",data["tag"]) )
            paramfile.write("NUMERIC_PROPERTY %s %f\n" % ("log_ki",data["value"]) )
            paramfile.write("STRING_PROPERTY %s %s\n" %("ligand_id",data["sdf_record"]))

def append_name_tag_to_params(data):
    if data["output_path"] != None:
        with open(data["output_path"],'a') as paramfile:
            paramfile.write("STRING_PROPERTY %s %s\n" %("ligand_id",data["sdf_record"]))    

def process_input_sdf(data):
    ligand_name = data["ligand_name"]
    input_sdf_path = data["filename"]
    output_dir = data["output_dir"]
    output_path = make_param_file(input_sdf_path,ligand_name,output_dir)
    data["output_path"] = output_path
    if "tag" in data and "value" in data:
        append_activity_tag_to_params(data)
    elif "output_path" in data:
        append_name_tag_to_params(data)
    return data
  
def set_up_crossed_ligands(active_list):
    new_inactive_list = []
    for active in active_list:
        for inactive in active_list:
            if active == inactive:
                continue
            
            inactive_record = inactive
            inactive_record["tag"] = active["tag"]
            inactive_record["value"] = 0.0
            new_inactive_list.append(inactive_record)
    return new_inactive_list
    
def add_ligand_names(ligand_path_data,rosetta_database):
    character_space = [
        'A','B','C','D','E',
        'F','G','H','I','J',
        'K','L','M','N','O',
        'P','Q','R','S','T',
        'U','V','W','X','Y',
        'Z','0','1','2','3',
        '4','5','6','7','8','9' ]
    name_size = required_ligand_name_length(len(ligand_path_data))
    names = itertools.product(character_space,repeat=name_size)
    disallowed_ligands =get_disallowed_ligands(rosetta_database)
    for index,record in enumerate(ligand_path_data):
        while True:
            ligand_name = names.next()
            ligand_name = "".join(ligand_name)
            if ligand_name not in disallowed_ligands:
                break
        
        ligand_path_data[index]["ligand_name"] = ligand_name
        
    return ligand_path_data


    
def init_options():
    usage = "%prog -jn ligands.db3 params_output_dir"
    parser=OptionParser(usage)
    parser.add_option("-j",dest="nprocs",default=2)
    parser.add_option("--database",dest="database",help="path to rosetta database",default="")
    parser.add_option("--path_to_params",dest="params_path",help="path to the molfile_to_params script", default="")
    parser.add_option("--make_all_ligands",dest="all_ligands",help="make params for all ligands even if no activities are present",default=False, action="store_true")
    parser.add_option("--cross_actives",dest="cross_actives",help="cross all active params as inactives",default=False,action="store_true")
    parser.add_option("--system_list",dest="system_list",help="only produce ligands for the specified systems",default=None)
    return parser
    
if __name__ == "__main__":
    
    (options,args) = init_options().parse_args()
    
    database_name = args[0]
    output_dir = args[1]
    path_to_params = options.params_path
    processor_pool = Pool(int(options.nprocs))
    
    setup_dir_for_filehash(output_dir)
    
    if options.system_list != None:
        with open(options.system_list) as system_file:
            system_list = []
            for line in system_file:
                line = line.rstrip()
                system_list.append(line)
    else:
        system_list = None
                
    active_information = []
    
    if options.all_ligands:
        activity_names = get_all_file_names(database_name,json_output=True)
    else:
        activity_names = get_file_names_with_activity_data(database_name,system_list=system_list)
    for record in activity_names:
        record["output_dir"] = output_dir
        active_information.append(record)
        
    if options.cross_actives:
        file_information = set_up_crossed_ligands(active_information)
    else:
        file_information = active_information
    add_ligand_names(file_information,options.database)
    
    complete_file_information = processor_pool.map(process_input_sdf,file_information)
    
    if not options.all_ligands:
        setup_params_schema(database_name)
        add_params_data(database_name,complete_file_information)
    
    
    
    
    