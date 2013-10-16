#!/usr/bin/env python

import sys
import glob
import subprocess

def run_bcl(input_file,output_dir):
    command = "bcl.exe model:ComputeStatistics "
    
    base_name = input_file.split("/")[-1].split(".")[0]
    flags = [
        "-input %s" % input_file,
        "-potency_cutoff 6.0",
        "-output_directory %s" % output_dir,
        "-table_name %s_stats.txt" %base_name
    ]
    
    command += " ".join(flags)
    
    subprocess.call(command,shell=True)
    

input_dir = sys.argv[1]
output_dir = sys.argv[2]

for path in glob.glob(input_dir+"/*"):
    run_bcl(path,output_dir)