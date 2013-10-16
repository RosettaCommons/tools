#!/usr/bin/env python

'''Given a set of molfiles scored by Rosetta, create a seperate BCL dataset for each molfile. 

Author: Sam DeLuca'''

import subprocess
import os
import sys
import glob
from optparse import OptionParser
from multiprocessing import Pool

def sort_sdf(job):
    bcl_path,input_path,output_path,delete_originals = job
    bcl_command = bcl_path
    
    flags = [
        "molecule:Reorder",
        "-input_filenames %s" % input_path,
        "-output %s" % output_path,
        "-sort 'MiscProperty(\"interface_delta_X\")'"
    ]
    
    bcl_command += " "+ " ".join(flags)
    
    subprocess.call(bcl_command, shell=True)
    
    if delete_originals:
        os.unlink(input_path)
    

def init_options():
    usage = "%prog --bcl_path=/path/to/bcl.exe input_directory/ output_directory/"
    parser=OptionParser(usage)
    parser.add_option("--bcl_path",dest="bcl_path",help="path to bcl.exe",default="bcl.exe")
    parser.add_option("-j",dest="nprocs",help="Number of cpus to use",default=1)
    parser.add_option("--delete_originals",dest="delete_originals",help="Delete the original input files", default=False, action="store_true")
    return parser
    
if __name__ == "__main__":
    
    options,args = init_options().parse_args()
    input_dir = args[0]
    output_dir = args[1]
    
    job_list = []
    for molfile_path in glob.glob(input_dir+"/*"):
        base_name = molfile_path.split("/")[-1].split(".")[0]
        output_path = "%s/%s.sdf.gz" % (output_dir,base_name)
        job_list.append( (options.bcl_path,molfile_path,output_path,options.delete_originals) )
        
    job_pool = Pool(int(options.nprocs))
    job_pool.map(sort_sdf,job_list)
    job_pool.close()
    #sort_sdf(options.bcl_path,molfile_path,output_path)
    
    