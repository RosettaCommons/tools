#!/usr/bin/env python

'''Given an sdf file scored by rosetta, score with the bcl and output a csv table

Author: Sam DeLuca'''

import subprocess
import os
import sys
import glob
from optparse import OptionParser
from random import shuffle

def score_sdf(bcl_path,input_filename,model_dir,output_table,nprocs,additional_score_list):
    
    additional_score_list.append("predicted_activity")
    formatted_score_terms = ["'Cached(%s)'" % term for term in additional_score_list]
    
    bcl_command = bcl_path
    
    flags = [
        "molecule:Properties",
        "-input_filenames %s" % input_filename,
        "-output_table %s" % output_path,
        "-add 'Mean(PredictedActivity(storage=File(directory=%s,prefix=model)))'" % model_dir,
        "-rename 'Mean(PredictedActivity(storage=File(directory=%s,prefix=model)))' predicted_activity" % model_dir,
        "-tabulate %s" % " ".join(formatted_score_terms),
    ]
    
    if nprocs > 1:
        flags.append("-scheduler PThread %d" % nprocs)
    
    bcl_command += " " + " ".join(flags)
    
    subprocess.call(bcl_command, shell=True)

def init_options():
    usage = "%prog --bcl_path=/path/to/bcl.exe --model_dir=ann_models input_file.sdf output_file.csv"
    parser=OptionParser(usage)
    parser.add_option("-j",dest="nprocs",help="number of cpus to use",default=1)
    parser.add_option("--bcl_path",dest="bcl_path",help="path to bcl.exe",default="bcl.exe")
    parser.add_option("--model_dir",dest="model_dir",help="directory with bcl models",default="")
    parser.add_option("--other_scores",dest="other_scores",help="additional scores to tabluate, comma seperated",default="")
    return parser
    
if __name__ == "__main__":
    
    options,args = init_options().parse_args()
    input_path = args[0]
    output_path = args[1]
    
    if options.other_scores != "":
        other_scores = options.other_scores.split(",")
    else:
        other_scores = []
    
    score_sdf(options.bcl_path,input_path,options.model_dir,output_path,int(options.nprocs),other_scores)
    