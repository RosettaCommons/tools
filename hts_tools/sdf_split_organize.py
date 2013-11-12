#!/usr/bin/env python
#definitely not the best script ive ever written

'''Split up an SDF file such that conformers of each individual ligand are in seperate files.
Output filenames will be hashed and split into directories by the first two digits of the hash
a path map with the original ligand ID and the new file path is output.

Author: Sam DeLuca'''

from optparse import OptionParser
import os
import re
import uuid 
import hashlib

def make_path(directory,name):
    sha1 = hashlib.sha1()
    sha1.update(name)
    digest = sha1.hexdigest()
    subdir = digest[0:2]
    path = "%(dir)s/%(subdir)s/%(digest)s.mol" % {"dir":directory,"subdir":subdir,"digest":digest}
    full_subdir = "%(dir)s/%(subdir)s" % {"dir":directory,"subdir":subdir}
    return (full_subdir,path)

def main():
    usage = "%prog infile outdir path_map"
    parser=OptionParser(usage)
    parser.add_option("-d",dest="descriptor",default="")
    parser.add_option("--make_uuid_names",dest="uuid",default=False, action="store_true")
    (options,args) = parser.parse_args()
    if len(args) != 3:
        parser.error("specify an input sdf file, an output directory, and a path to a path map file")
    infile = open(args[0],'r')
    path_map = open(args[2],'w')
    line_counter=1
    current_record = ""
    outfile = None
    record_count = 0
    line_buffer = []
    path_tuples = set()
    found_descriptor = False
    
    #write header
    path_map.write("ligand_id,filename\n")
    path_map.write("string,string\n")
    
    for line in infile:
        subdir,output_path = make_path(args[1],current_record)
        if not os.path.exists(subdir):
            os.mkdir(subdir)
        if(found_descriptor):
            line=line.strip()
            line_buffer.append(line)
            current_record = line
            current_record = current_record.split(",")[-1].strip()
            print current_record
            found_descriptor = False
        if(line[0] is "$"):
            line_counter = 0
            record_count += 1
            if(os.path.exists(output_path)):
                outfile = open(output_path,'a')
            else:
                outfile = open(output_path,'w')
            outfile.writelines(line_buffer)
            if(outfile):
                outfile.write("$$$$\n")
                outfile.close()
                path_tuples.add( (current_record,output_path) )
                #path_map.write(current_record+","+output_path+"\n")
            line_buffer = []
        elif(line_counter is 1 and options.descriptor is ""):
            if options.uuid:
                current_record = uuid.uuid1().hex
            else:
                current_record = line.rstrip()
                current_record = current_record.split(",")[-1].strip()
            if(os.path.exists(output_path)):
                outfile = open(output_path,'a')
            else:
                outfile = open(output_path,'w')
            #outfile.write(line)
            if options.uuid:
                line_buffer.append(current_record+"\n")
            else:
                line_buffer.append(line)
        elif(line_counter is not 1 and options.descriptor is not "" and line[0] is ">"):
            if(re.search(options.descriptor,line)):
                found_descriptor=True
                line_counter+=1
                line_buffer.append(line)
                continue
        else:
            line_buffer.append(line)
            #outfile.write(line)
        line_counter +=1
    print str(record_count) + " records processed"
    for current_record,output_path in path_tuples:
        path_map.write(current_record+","+output_path+"\n")
    path_map.close()

if __name__ == "__main__":
    main()
