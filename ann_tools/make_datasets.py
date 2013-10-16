#!/usr/bin/env python

'''Given a set of molfiles scored by Rosetta, create a seperate BCL dataset for each molfile. 

Author: Sam DeLuca'''

import subprocess
import os
import sys
import glob
from optparse import OptionParser
from random import shuffle
from multiprocessing import Pool
def make_data(job):
    
    bcl_path,chunk_files,output_path,input_features,output_features,first_n = job
    bcl_command = bcl_path
    
    flags = [
        "descriptor:GenerateDataset",
        "-feature_labels %s" % input_features,
        "-result_labels %s" % output_features,
        "-output %s" % output_path,
        "-scheduler PThread 8",
    ]
    
    source = "-source 'Balanced("
    for path in chunk_files:
        if first_n == 0:
            source+= "SdfFile(filename=%s)," % path
        else:
            source+= "Chunk(chunks=\"[0,%d)\",dataset=SdfFile(filename=%s))," %(first_n,path)
    source = source.strip(",") #get rid of the last trailing ,
    source += ")'"
    flags.append(source)
    bcl_command += " " + " ".join(flags)
    
    subprocess.call(bcl_command, shell=True)
    
    
def chunks(list, count):
    return [list[i:i+count] for i in range(0, len(list), count)]
    
def init_options():
    usage = "%prog --bcl_path=/path/to/bcl.exe path_list.txt output_directory"
    parser=OptionParser(usage)
    parser.add_option("-j",dest="nprocs",help="number of cpus to use",default=1)
    parser.add_option("--bcl_path",dest="bcl_path",help="path to bcl.exe",default="bcl.exe")
    parser.add_option("--input_features",dest="input_features",help="input feature label file",default="")
    parser.add_option("--output_features",help="output feature label file",default="")
    parser.add_option("--first_n",dest="first_n",help="only include the first n compounds from each input file",default=0)
    parser.add_option("--chunk_count",dest="n_chunks",help="process input sdf files into n chunks",default=1)
    parser.add_option("--output_prefix",dest="prefix",help="prefix for output files",default="dataset")
    return parser
    
if __name__ == "__main__":
    
    options,args = init_options().parse_args()
    input_list = args[0]
    output_dir = args[1]
    
    with open(input_list) as infile:
        molfile_list = [molfile.strip() for molfile in infile]
    shuffle(molfile_list)
    chunk_list = chunks(molfile_list,len(molfile_list)/int(options.n_chunks))
    
    chunk_counter = 1
    job_list = []
    for chunk in chunk_list:
        output_path = "%s/%s_%03d.bin" % (output_dir,options.prefix,chunk_counter)
        chunk_counter += 1
        job_list.append( (options.bcl_path,chunk,output_path,options.input_features,options.output_features,int(options.first_n)) )
    job_pool = Pool(int(options.nprocs))
    job_pool.map(make_data,job_list)
    job_pool.close()
    
    