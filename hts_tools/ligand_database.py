import sqlite3


def setup_input_schema(db_name,header):
    '''given a header vector and a db_name, set up a basic schema'''
    
    schema_fields = ",".join([key + " "+ header[key] for key in header])
    
    schema_string = "CREATE TABLE IF NOT EXISTS sdf_input_data (record_id INTEGER PRIMARY KEY AUTOINCREMENT," + schema_fields + ")"
    
    connection = sqlite3.connect(db_name)
    cursor = connection.cursor()
    cursor.execute(schema_string)
    connection.commit()
    connection.close()
    
def setup_tag_schema(db_name):
    '''given a db name, set up the tag schema'''
    
    schema_string = "CREATE TABLE IF NOT EXISTS activity_tags (tag_id INTEGER PRIMARY KEY AUTOINCREMENT,record_id INTEGER REFERENCES sdf_input_data(record_id), tag string,value float)"
    connection = sqlite3.connect(db_name)
    cursor = connection.cursor()
    cursor.execute(schema_string)
    connection.commit()
    connection.close()
    
def write_data(db_name,table,columns,data_list):
    '''Given a data map, list of columns and a table name, write data to the database'''
    column_inserts = "(" + ",".join(["?" for x in range(len(columns))]) + ")"
    column_names = "(" + ",".join(columns) + ")"
    insert_string = "INSERT INTO " + table + " " + column_names + " VALUES " + column_inserts
    connection = sqlite3.connect(db_name)
    cursor = connection.cursor()
    for record in data_list:
        data_list = [record[key] for key in columns]
        cursor.execute(insert_string,data_list)
    connection.commit()
    connection.close()
    
def write_tag_data(db_name,data_list):
    insert_string = "INSERT INTO activity_tags (record_id,tag,value) VALUES (?,?,?)"
    query_string = "SELECT record_id FROM sdf_input_data WHERE ligand_id = ?"
    
    connection = sqlite3.connect(db_name)
    cursor = connection.cursor()
    
    query_record_id_map = {}
    #print "getting ids for",len(data_list),"values"
    #for record in data_list:
    #    ligand_id = record["ligand_id"]
    #    for row in cursor.execute(query_string,(ligand_id,)):
    #        query_record_id_map[ligand_id] =row[0]
    #print "writing activity tags"
    skip_count = 0
    for record in data_list:
        ligand_id = record["ligand_id"]
        record_id = None
        cursor.execute(query_string,(ligand_id,))
        try:
            record_id = cursor.fetchone()[0]
        except TypeError:
            print "skipping",ligand_id
            skip_count += 1
            continue
        print "found",ligand_id
        tag = record["tag"]
        value = record["value"]
        cursor.execute(insert_string,(record_id,tag,value))
        connection.commit()
    print "skippped",skip_count
    connection.close()
    

def get_all_file_names(db_name):
    select_string = "SELECT record_id,filename FROM sdf_input_data"
    connection = sqlite3.connect(db_name)
    cursor = connection.cursor()
    for record_id,filename in cursor.execute(select_string):
        yield (record_id,filename)
        
    connection.close()
    
def update_filenames(db_name, filename_update_data):
    update_string = "UPDATE sdf_input_data SET filename = ? WHERE record_id = ?"
    connection = sqlite3.connect(db_name)
    cursor = connection.cursor()
    for record_id,path in filename_update_data:
        cursor.execute(update_string,(path,record_id))
        connection.commit()
    connection.close()