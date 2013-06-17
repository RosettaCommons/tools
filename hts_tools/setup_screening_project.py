#!/usr/bin/env python

'''This is script will construct the initial database used for managing SDF and param files
for rosettaHTS screening projects.

Author: Sam DeLuca'''

from optparse import OptionParser
import csv
import sys
from ligand_database import *
from os.path import exists

def parse_input_file(input_path):
    '''Parse the CSV input file into a list of dicts.  Return the header and '''
    header = []
    header_schema = {}
    data_list = []
    with open(input_path,'r') as csvfile:
        for index, row in enumerate(csv.reader(csvfile)):
            if index == 0:
                header = row
                if "filename" not in header:
                    sys.exit("ERROR: input file header must have a column labeled 'filename'")
                if "ligand_id" not in header:
                    sys.exit("ERROR: input file header must have a column labeled 'ligand_id'")
            elif index == 1:
                for column,datatype in zip(header,row):
                    header_schema[column] = datatype
            else:
                row_map = {}
                for column, value in zip(header,row):
                    try:
                        row_map[column] = float(value)
                    except ValueError:
                        row_map[column] = value
                data_list.append(row_map)
        return (header_schema, data_list)
        
def all_files_exist(data_list):
    for record in data_list:
        if not exists(record["filename"]):
            return (False, record)
    return (True, None)

def init_options():
    usage = "%prog input_file.csv output.db3"
    parser=OptionParser(usage)
    
    return parser

if __name__ == "__main__":
    
    options,args = init_options().parse_args()
    if len(args) != 2:
        parser.error("you must specify both an input csv file and an output database")
    input_file = args[0]
    output_db = args[1]
    
    #parse input file
    header, data_list = parse_input_file(input_file)
    
    status, record = all_files_exist(data_list)
    if status == False:
        sys.exit("ERROR: %(filename)s does not exist" % {"filename" : record["filename"]})
    
    setup_input_schema(output_db,header)
    
    write_data(output_db,"sdf_input_data",header.keys(),data_list)


    