#!/usr/bin/env python

'''This is script will construct the initial database used for managing SDF and param files
for rosettaHTS screening projects.

Author: Sam DeLuca'''

from optparse import OptionParser
import csv
import sys

def parse_input_file(input_path):
    '''Parse the CSV input file into a list of dicts.  Return the header and '''
    header = []
    data_list = []
    with open(input_path,'r') as csvfile:
        for index, row in enumerate(csv.reader(csvfile)):
            if index == 0:
                header = row
                if "filename" not in header:
                    sys.exit("ERROR: input file header must have a column labeled 'filename'")
                if "ligand_id" not in header:
                    sys.exit("ERROR: input file header must have a column labeled 'ligand_id'")
                    
            else:
                row_map = {}
                for header, value in zip(header,row):
                    try:
                        row_map[header] = float(value)
                    except ValueError:
                        row_map[header] = value
                data_list.append(row_map)
        return (header, data_list)
        
def 


def init_options():
    usage = "%prog input_file.csv"
    parser=OptionParser(usage)
    
    return parser

if __name__ == "__main__":
    
    options,args = parser.parse_args()
    input_file = args[1]
    
    #parse input file
    header, data_list = parse_input_file(input_file)
    
    #validate paths
    
    #create schema
    
    #insert data into schema
    