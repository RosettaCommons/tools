import csv
import sys

def parse_input_file(input_path,required_headers):
    '''Parse the CSV input file into a list of dicts.  Return the header and '''
    header = []
    header_schema = {}
    data_list = []
    with open(input_path,'r') as csvfile:
        for index, row in enumerate(csv.reader(csvfile)):
            if index == 0:
                header = row
                for req in required_headers:
                    if req not in header:
                        sys.exit("ERROR: input file header must have a column labeled "+req)
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