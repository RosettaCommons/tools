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
    print column_names,column_inserts
    insert_string = "INSERT INTO " + table + " " + column_names + " VALUES " + column_inserts
    print insert_string
    connection = sqlite3.connect(db_name)
    cursor = connection.cursor()
    for record in data_list:
        data_list = [record[key] for key in columns]
        cursor.execute(insert_string,data_list)
    connection.commit()
    connection.close()