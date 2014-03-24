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
    
def setup_metadata_schema(db_name):
    '''given a db name, set up the metadata schema'''
    
    schema_string = "CREATE TABLE IF NOT EXISTS metadata (tag_id INTEGER PRIMARY KEY AUTOINCREMENT,record_id INTEGER REFERENCES sdf_input_data(record_id), tag string,value float)"
    connection = sqlite3.connect(db_name)
    cursor = connection.cursor()
    cursor.execute(schema_string)
    connection.commit()
    connection.close()
    
def setup_params_schema(db_name):
    '''given a db name, setup the params schema'''
    schema_string = "CREATE TABLE IF NOT EXISTS params_data (param_id INTEGER PRIMARY KEY AUTOINCREMENT, tag_id INTEGER REFERENCES activity_tags(tag_id),ligand_name string,filename string)"
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
    
def write_tag_data(db_name,data_list,mode="tag"):
    '''given a deata list, write activity tag information'''
    if mode =="tag":
        insert_string = "INSERT INTO activity_tags (record_id,tag,value) VALUES (?,?,?)"
    elif mode== "metadata":
        insert_string = "INSERT INTO metadata (record_id,tag,value) VALUES (?,?,?)"
    query_string = "SELECT ligand_id,record_id FROM sdf_input_data"
    
    
    
    connection = sqlite3.connect(db_name)
    cursor = connection.cursor()
    query_record_id_map = {}
    
    for ligand_id,record_id in cursor.execute(query_string):
        query_record_id_map[ligand_id] = record_id
    
    skip_count = 0
    
    #we're about to insert a lot of rows, turn some knobs to help
    cursor.execute("PRAGMA synchronous=OFF")
    cursor.execute("PRAGMA journal_mode=MEMORY")
    
    for record in data_list:
        ligand_id = record["ligand_id"]
        try:
            record_id = query_record_id_map[ligand_id]
        except KeyError:
            skip_count += 1
            #print "skipping",ligand_id
            continue
        #print "found",ligand_id
        tag = record["tag"]
        value = record["value"]
        cursor.execute(insert_string,(record_id,tag,value))
    connection.commit()
    print "skippped",skip_count
    
    #go back to safer settings
    cursor.execute("PRAGMA synchronous=ON")
    cursor.execute("PRAGMA journal_mode=DELETE")
    
    connection.close()
    

def get_all_file_names(db_name,only_tagged=False,json_output=False):
    '''Generator producing file names for every input sdf'''
    if only_tagged:
        select_string = "SELECT sdf_input_data.ligand_id,filename FROM sdf_input_data JOIN activity_tags ON sdf_input_data.record_id == activity_tags.record_id"
    else:
        select_string = "SELECT ligand_id,filename FROM sdf_input_data"
    connection = sqlite3.connect(db_name)
    cursor = connection.cursor()
    
    if json_output:
        for record_id, filename in cursor.execute(select_string):
            yield {
                "sdf_record" : record_id,
                "filename" : filename
            }
    else:
        for record_id,filename in cursor.execute(select_string):
            yield (record_id,filename)
        
    connection.close()
    
def get_file_names_with_activity_data(db_name,system_list = None):
    '''generator producing sdf filenames and activity tag data'''
    select_string = "SELECT sdf_input_data.ligand_id,activity_tags.tag_id,filename,tag,value FROM sdf_input_data JOIN activity_tags ON sdf_input_data.record_id == activity_tags.record_id"
    connection = sqlite3.connect(db_name)
    cursor = connection.cursor()
    for sdf_record,tag_id,filename,tag,value in cursor.execute(select_string):
        if system_list != None and tag in system_list:
            yield {
                "sdf_record" : sdf_record,
                "tag_id" : tag_id,
                "filename" : filename,
                "tag" : tag,
                "value" : value
            }
        elif system_list == None:
            yield {
                "sdf_record" : sdf_record,
                "tag_id" : tag_id,
                "filename" : filename,
                "tag" : tag,
                "value" : value
            }
    connection.close()
    
def get_params_information(db_name):
    '''Generator produce params file path information'''
    select_string = "SELECT tag, filename FROM params_data JOIN activity_tags ON params_data.tag_id = activity_tags.tag_id"
    connection = sqlite3.connect(db_name)
    cursor = connection.cursor()
    
    for tag,filename in cursor.execute(select_string):
        yield {
            "tag" : tag,
            "params_file" : filename
        }
    connection.close()
        
def get_params_sdf_mapping(db_name):
    '''Generator to produce mapping between params name and sdf path'''
    select_string = "SELECT ligand_name,sdf_input_data.filename FROM params_data JOIN activity_tags ON params_data.tag_id = activity_tags.tag_id JOIN sdf_input_data ON activity_tags.record_id = sdf_input_data.record_id"
    
    connection = sqlite3.connect(db_name)
    cursor = connection.cursor()
    for  name, filename in cursor.execute(select_string):
        yield {
            "ligand_name" : name,
            "file_name" : filename
        }
    connection.close()

def get_sdf_path_from_ligand_name(db_name,ligand_name):
    '''Given a ligand_name, return the filename of the sdf file'''
    select_string = "SELECT sdf_input_data.filename FROM params_data JOIN activity_tags ON params_data.tag_id = activity_tags.tag_id JOIN sdf_input_data ON activity_tags.record_id = sdf_input_data.record_id WHERE ligand_name = ?"
    connection = sqlite3.connect(db_name)
    cursor = connection.cursor()
    cursor.execute(select_string,(ligand_name,))
    filename = cursor.fetchone()[0]
    connection.close()
    return filename

def add_params_data(db_name,data_list):
    '''add information about a new params filename'''
    connection = sqlite3.connect(db_name)
    cursor = connection.cursor()
    insert_string = "INSERT INTO params_data  (tag_id,ligand_name,filename) VALUES (?,?,?)"
    for record in data_list:
        tag_id = record["tag_id"]
        ligand_name = record["ligand_name"]
        filename = record["output_path"]
        if filename == None:
            continue
        cursor.execute(insert_string,(tag_id,ligand_name,filename))
        connection.commit()
    connection.close()
        
    
def update_filenames(db_name, filename_update_data):
    '''update an sdf file path'''
    update_string = "UPDATE sdf_input_data SET filename = ? WHERE record_id = ?"
    connection = sqlite3.connect(db_name)
    cursor = connection.cursor()
    for record_id,path in filename_update_data:
        cursor.execute(update_string,(path,record_id))
        connection.commit()
    connection.close()