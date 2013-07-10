#!/usr/bin/env python

'''
Extract data from the screening_features reporter into a json file. For now we're just supporting mysql, but this script uses sqlalchemy and can thus be easily extended to support sqlite3, postgres, etc
Author: Sam DeLuca'''

from __future__ import with_statement
import getpass
import sys
import sqlalchemy
from sqlalchemy import and_
import sqlalchemy.orm
from optparse import OptionParser
try:
    import json
except ImportError:
    import simplejson as json

def page_query(q):
	'''Break up huge queries into chunks to avoid insane memory usage'''
	chunk = 100000
	offset = 0
	while True:
		r = False
		for elem in q.limit(chunk).offset(offset):
			r = True
			yield elem
		offset += chunk
		if not r:
			break
			
def make_table(tablename,engine):
	'''Helper function for making table objects via introspection'''
	meta = sqlalchemy.MetaData(engine)
	return sqlalchemy.Table(tablename,meta,autoload=True,autoload_with=engine)


def get_ligand_data(SessionMaker,engine,batch_descriptor):
    "get data records for each ligand"
    query_session = SessionMaker()
    screening_features = make_table('screening_features',engine)
    structures = make_table('structures',engine)
    batches = make_table('batches',engine)
    
    ligand_data_list = []
    struct_id_query = query_session.query(
        structures.columns.struct_id).\
        join(batches,structures.columns.batch_id == batches.columns.batch_id).\
        filter(batches.columns.description == batch_descriptor)
        
    struct_ids = set()
    print "Getting struct_ids"
    for struct_id in page_query(struct_id_query):
        struct_ids.add(struct_id[0])
    print "Got",len(struct_ids),"valid struct_ids"
    for struct_id in struct_ids:
        data_query = query_session.query(
                screening_features.columns.name3,
                screening_features.columns.group_name,
                screening_features.columns.descriptor_data).\
                filter(screening_features.columns.struct_id == struct_id).one()
        name3,group_name,descriptor_data_string = data_query
        data_entry = {}
        data_entry["struct_id"] = struct_id
        data_entry["name3"] = name3
        data_entry["group_name"] = group_name
        try:
            descriptor_data = json.loads(descriptor_data_string)
        except json.decoder.JSONDecodeError:
            print "Can't decode data for struct_id",struct_id
            print descriptor_data_string
            continue
        for key in descriptor_data:
            data_entry[key] = descriptor_data[key]
        #ligand_data_list.append(data_entry)
        yield data_entry
    #return ligand_data_list
    query_session.close()
    
def init_options():
    usage = "%prog  --batch-description=batch_desc  output_dir/"
    parser = OptionParser(usage)
    parser.add_option("--batch-description",dest="batch_description",help="select the batch description tag to process data from", default="")
    parser.add_option("--host",dest="host",help="Server host",default="")
    parser.add_option("--username",dest="username",help="Database user name",default="")
    parser.add_option("--database-name",dest="db_name",help="The name of the database to connect to",default="")
    parser.add_option("--chunk-size",dest="chunk_size",help="Number of records per output_file (default 10000)",default="10000")
    return parser
    
if __name__ == "__main__":
    
    parser = init_options()
    (options,args) = parser.parse_args()
    
    output_directory = args[0]
    
    if options.batch_description == "":
        parser.exit("You must specify --batch-description")
    if options.host == "":
        parser.exit("You must specify --host")
    if options.username == "":
        parser.exit("You must specify --username")
    if options.db_name == "":
        parser.exit("You must specify --database-name")
    if len(args) != 1:
        parser.exit("You must specify an output filename")
    
    password = getpass.getpass("Enter Password for user " + options.username+ ": ")
    
    mysql_engine_string = "mysql+mysqldb://%s:%s@%s/%s"
   
    
    engine_data = (
        options.username,
        password,
        options.host,
        options.db_name
    )
    
    engine = sqlalchemy.create_engine(mysql_engine_string % engine_data,pool_recycle=3600)
    SessionMaker = sqlalchemy.orm.sessionmaker(bind=engine)
    chunk_size = int(options.chunk_size)
    chunk_number = 1
    
    partial_data = []
    for ligand_data in get_ligand_data(SessionMaker,engine,options.batch_description):
        partial_data.append(ligand_data)
        if len(partial_data) == chunk_size:
            output_filename = "%s/%s_%04d.js" % (output_directory,options.batch_description,chunk_number)
            with open(output_filename,'w') as output:
                json.dump(partial_data,output)
            partial_data = []
            chunk_number += 1
    