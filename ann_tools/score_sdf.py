#!/usr/bin/env python

'''Given an sdf file scored by rosetta, score with the bcl and output a csv table

Author: Sam DeLuca'''

import subprocess
import os
import sys
import glob
from optparse import OptionParser
from multiprocessing import Pool


def score_sdf(job):
    bcl_path,input_filenames,model_dir,output_path,additional_score_list = job
    additional_score_list.append("predicted_activity")
    formatted_score_terms = ["'Cached(%s)'" % term for term in additional_score_list]
    formatted_score_terms.append("'Cached(Weight)'")
    bcl_command = bcl_path
    
    flags = [
        "molecule:Properties",
        "-input_filenames %s" % " ".join(input_filenames),
        "-output_table %s" % output_path,
        "-add 'Weight' 'Mean(PredictedActivity(storage=File(directory=%s,prefix=model)))'" % model_dir,
        "-rename 'Mean(PredictedActivity(storage=File(directory=%s,prefix=model)))' predicted_activity" % model_dir,
        "-tabulate %s" % " ".join(formatted_score_terms),
    ]
    
    bcl_command += " " + " ".join(flags)
    
    subprocess.call(bcl_command, shell=True)
    
def chunks(l, n):
    for i in xrange(0, len(l), n):
        yield l[i:i+n]
        
def init_options():
    usage = "%prog --bcl_path=/path/to/bcl.exe --model_dir=ann_models input_list.txt output_prefix"
    parser=OptionParser(usage)
    parser.add_option("-j",dest="nprocs",help="number of cpus to use",default=1)
    parser.add_option("--bcl_path",dest="bcl_path",help="path to bcl.exe",default="bcl.exe")
    parser.add_option("--model_dir",dest="model_dir",help="directory with bcl models",default="")
    parser.add_option("--other_scores",dest="other_scores",help="additional scores to tabluate, comma seperated",default="")
    return parser
    
if __name__ == "__main__":
    
    options,args = init_options().parse_args()
    input_list = args[0]
    output_prefix = args[1]
    
    if options.other_scores != "":
        other_scores = options.other_scores.split(",")
    else:
        other_scores = []
    
    input_paths = []
    with open(input_list) as infile:
        for path in infile:
            input_paths.append(path.strip())
    
    jobs = []
    for index,chunk in enumerate(chunks(input_paths,len(input_paths)/int(options.nprocs))):
        jobs.append( (options.bcl_path,chunk,options.model_dir,output_prefix+"_"+str(index)+".csv",other_scores) )
    
    workers = Pool(int(options.nprocs))
    workers.map(score_sdf,jobs)
    workers.close()
    workers.join()
    #score_sdf(options.bcl_path,input_path,options.model_dir,output_path,other_scores)
    