#!/usr/bin/env python

'''Given a database created with setup_screening_project, add hydrogens to every ligand file and compress. You need to have corina in your path to use this script.   If you don't have corina, modify add_hydrogens to use your hydrogen generation method

Author: Sam DeLuca'''

import tempfile
import subprocess
import os
import sys
from ligand_database import *
from optparse import OptionParser
from multiprocessing import Pool
import shutil 

def already_processed(input_file):
    '''Figure out if corina has already processed a file'''
    in_handle = open(input_file)
    for i, line in enumerate(in_handle):
        if i == 1:
            fields = line.split()
            if fields[0] == "SDcorina":
                #print "skipping", input_file.split("/")[-1]
                in_handle.close()
                return True
                
            else:
                print "processing",input_file.split("/")[-1]
                in_handle.close()
                return False
                
def add_hydrogens(input_file):
    '''Run corina to add hydrogens'''
    if already_processed(input_file):
        return
    out_file = input_file+".tmp"
    corina_command = "corina -d wh,no3d %(infile)s %(outfile)s"
    #print corina_command % {"infile": input_file,"outfile" : out_file}
    subprocess.call(corina_command % {"infile": input_file,"outfile" : out_file}, shell=True)
    os.rename(out_file,input_file)
    
def init_options():
    usage = "%prog -jn ligands.db3"
    parser=OptionParser(usage)
    parser.add_option("-j",dest="nprocs",default=2)
    parser.add_option("--only_tagged",dest="only_tagged",help="Only add hydrogens for ligands with activity tags",default=False,action="store_true")
    return parser
    
if __name__ == "__main__":
    options,args = init_options().parse_args()
    if len(args) != 1:
        parser.error("you must specify both an input database")
    database_path = args[0]
    
    if not os.path.exists(database_path):
        sys.exit(database_path+" does not exist")

    processor_pool = Pool(int(options.nprocs))

    filenames = []
    for data in get_all_file_names(database_path,only_tagged=options.only_tagged):
        if data[1].split(".")[-1] == "gz":
            sys.exit("Some or all of the files in the database are compressed.  Corina doesn't work on gzipped files")
        filenames.append(data[1])
    
    new_filenames = processor_pool.map(add_hydrogens,filenames)