#!/usr/bin/env python

'''Given a database created with setup_screening_project, gzip compress every ligand file and update the paths stored in the database to match.

Author: Sam DeLuca'''

import sys
import gzip
import os
from ligand_database import *
from optparse import OptionParser
from multiprocessing import Pool

def compress_file(path_to_compress):
    zip_path = path_to_compress+".gz"
    with open(path_to_compress,'r') as infile:
        outfile = gzip.open(zip_path,'wb')
        outfile.writelines(infile)
        outfile.close()
    os.unlink(path_to_compress)
    return zip_path
    
def process_file(data):
    record_id, input_path = data
    output_path = compress_file(input_path)
    return (record_id,output_path)

def init_options():
    usage = "%prog -jn ligands.db3"
    parser=OptionParser(usage)
    parser.add_option("-j",dest="nprocs",default=2)
    return parser

if __name__ == "__main__":
    options,args = init_options().parse_args()
    if len(args) != 1:
        parser.error("you must specify both an input database")
    database_path = args[0]
    if not os.path.exists(database_path):
        sys.exit(database_path+" does not exist")
    
    processor_pool = Pool(int(options.nprocs))
    
    file_records = []
    for data in get_all_file_names(database_path):
        file_records.append(data)
        
    new_filenames = processor_pool.map(process_file,file_records)
    
    update_filenames(database_path,new_filenames)