#!/usr/bin/env python

'''This is script will construct the initial database used for managing SDF and param files
for rosettaHTS screening projects.

Author: Sam DeLuca'''

from optparse import OptionParser
import csv
import sys
from ligand_database import *
from os.path import exists
from hts_util import parse_input_file
        
def all_files_exist(data_list):
    for record in data_list:
        if not exists(record["filename"]):
            return (False, record)
    return (True, None)

def init_options():
    usage = "%prog input_file.csv output.db3"
    
    parser=OptionParser(usage)
    parser.add_option("--no_verify",dest="no_verify",help="Don't check that paths in input csv exist before writing db",default=False, action="store_true")
    return parser

if __name__ == "__main__":
    
    options,args = init_options().parse_args()
    if len(args) != 2:
        parser.error("you must specify both an input csv file and an output database")
    input_file = args[0]
    output_db = args[1]
    
    #parse input file
    header, data_list = parse_input_file(input_file,["filename","ligand_id"])
    
    if not options.no_verify:
        status, record = all_files_exist(data_list)
        if status == False:
            sys.exit("ERROR: %(filename)s does not exist" % {"filename" : record["filename"]})
    
    setup_input_schema(output_db,header)
    
    write_data(output_db,"sdf_input_data",header.keys(),data_list)


    