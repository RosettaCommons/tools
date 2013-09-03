#!/usr/bin/env python

'''Given a set of molfiles scored by Rosetta, create a seperate BCL dataset for each molfile. 

Author: Sam DeLuca'''

import subprocess
import os
import sys
import glob
from optparse import OptionParser


def sort_sdf(bcl_path,input_path,output_path):
    bcl_command = bcl_path
    
    flags = [
        "molecule:Reorder",
        "-input_filenames %s" % input_path,
        "-output %s" % output_path,
        "-sort 'MiscProperty(\"interface_delta_X\")'"
    ]
    
    bcl_command += " "+ " ".join(flags)
    
    subprocess.call(bcl_command, shell=True)
    

def init_options():
    usage = "%prog --bcl_path=/path/to/bcl.exe input_directory/ output_directory/"
    parser=OptionParser(usage)
    parser.add_option("--bcl_path",dest="bcl_path",help="path to bcl.exe",default="bcl.exe")
    return parser
    
if __name__ == "__main__":
    
    options,args = init_options().parse_args()
    input_dir = args[0]
    output_dir = args[1]
    
    for molfile_path in glob.glob(input_dir+"/*"):
        base_name = molfile_path.split("/")[-1].split(".")[0]
        output_path = "%s/%s.sdf" % (output_dir,base_name)
        sort_sdf(options.bcl_path,molfile_path,output_path)
    
    