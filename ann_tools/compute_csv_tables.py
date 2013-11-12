#!/usr/bin/env python2.7

import subprocess
import os
import sys
import glob
import numpy
from optparse import OptionParser
from random import shuffle
from multiprocessing import Pool

def make_data(job):
    
    try:
        bcl_path,input_file,output_path,input_features,output_features,first_n = job
        bcl_command = bcl_path
    
        flags = [
            "descriptor:GenerateDataset",
            "-feature_labels %s" % input_features,
            "-result_labels %s" % output_features,
            "-output %s" % output_path,
        ]
    
        source = "-source 'Chunk(chunks=\"[0,%d)\",dataset=SdfFile(filename=%s))'" %(first_n,input_file)
        source = source.strip(",") #get rid of the last trailing ,
        flags.append(source)
        bcl_command += " " + " ".join(flags)
    
        if subprocess.call(bcl_command, shell=True,stdout=open(os.devnull,'wb')) != 0:
            print "Skipping",input_file,"due to some error"
            return None
    
        array_data = []
        with open(output_path) as outfile:
            for line in outfile:
                values = [float(x) for x in line.split(",")]
                array_data.append(values)
    
        means = numpy.mean(array_data,axis=0)
        stdevs = numpy.std(array_data[0:-1],axis=0)
    
        means[numpy.abs(means) < 1e-10] = 0
        stdevs[numpy.abs(stdevs) < 1e-10] = 0
    
        final_row = numpy.append(stdevs,means)
        os.unlink(output_path)
        return final_row
    except:
        print "there was some error in the process, bailing out"
        return None
    
            

def init_options():
    usage = "%prog --bcl_path=/path/to/bcl.exe path_list.txt output_directory"
    parser=OptionParser(usage)
    parser.add_option("-j",dest="nprocs",help="number of cpus to use",default=1)
    parser.add_option("--bcl_path",dest="bcl_path",help="path to bcl.exe",default="bcl.exe")
    parser.add_option("--input_features",dest="input_features",help="input feature label file",default="")
    parser.add_option("--output_features",help="output feature label file",default="")
    parser.add_option("--first_n",dest="first_n",help="only include the first n compounds from each input file",default=0)
    parser.add_option("--output_prefix",dest="prefix",help="prefix for output files",default="dataset")
    return parser

if __name__ == "__main__":
    
    options,args = init_options().parse_args()
    input_list = args[0]
    output_filename = args[1]
    
    with open(input_list) as infile:
        molfile_list = [molfile.strip() for molfile in infile]
    
    job_list = []
    for path in molfile_list:
        path = path.strip()
        base_name = path.split("/")[-1].split(".")[0]
        output_path = "%s_%s.csv" % (options.prefix,base_name)
        job_list.append( (options.bcl_path,path,output_path,options.input_features,options.output_features,int(options.first_n)) )
    
    job_pool = Pool(int(options.nprocs))
    processed_rows = job_pool.map(make_data,job_list)
    job_pool.close()
    
    with open(output_filename,'w') as outfile:
        for row in processed_rows:
            if row == None:
                continue
            formatted_line = ",".join([str(x) for x in row])
            outfile.write("%s\n"%formatted_line)
            