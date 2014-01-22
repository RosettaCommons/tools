#!/usr/bin/env python

'''Given a database created with setup_screening_project, add tags and activity data.  Each ligand can have multiple tags.
This script is used to add data about ligand activity to ligand records. 

Author: Sam DeLuca'''

import tempfile
import sys
from ligand_database import *
from optparse import OptionParser
from hts_util import parse_input_file

def init_options():
    usage = "%prog ligands.db3 tag file"
    parser=OptionParser(usage)
    parser.add_option("--mode",dest="mode",help="The type of tag to insert, should be 'metadata' or 'tag'. (Default: 'tag')",default="tag")
    return parser
    
if __name__ == "__main__":
    options,args = init_options().parse_args()
    if len(args) != 2:
        parser.error("you must specify both an input database and a tag file")
    database_path = args[0]
    tag_file = args[1]
    tag_file_schema,tag_file_list = parse_input_file(tag_file,["ligand_id","tag","value"])
    
    if options.mode == "tag":
        setup_tag_schema(database_path)
        write_tag_data(database_path,tag_file_list)
    elif option.mode == "metadata":
        setup_metadata_schema(database_path)
        write_tag_data(database_path,tag_file_list,mode="metadata")