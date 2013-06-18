import sqlite3


def setup_input_schema(db_name,header):
    '''given a header vector and a db_name, set up a basic schema'''
    
    schema_fields = ",".join([key + " "+ header[key] for key in header])
    
    schema_string = "CREATE TABLE sdf_input_data (record_id INTEGER PRIMARY KEY AUTOINCREMENT," + schema_fields + ")"
    
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