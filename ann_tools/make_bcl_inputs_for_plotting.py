#!/usr/bin/env python

import pandas as pd
import sys
import glob

def make_bcl_file(file_path,data):
    best_experimental = []
    best_predicted = []
    for name,group in data.groupby("Cached(ligand_id)"):
        top_row = group.sort(columns="Cached(total_score_X)")[0:10].mean(axis=0)
        experimental = top_row["Cached(log_ki)"]
        predicted = top_row["Cached(predicted_activity)"]
        if experimental == 1.0:
            experimental = 7.0
        best_experimental.append(experimental)
        best_predicted.append(predicted)
    
    record_count = len(best_experimental)
    
    with open(file_path,'w') as outfile:
        outfile.write("bcl::linal::Matrix<float>\n")
        outfile.write("%d\t2\n" % record_count)
        for experimental,predicted in zip(best_experimental,best_predicted):
            outfile.write( "%f\t%f\n" % (experimental,predicted) )

input_dir = sys.argv[1]
output_dir = sys.argv[2]

individual_frames = []
for input_path in glob.glob(input_dir+"/*.csv"):
    individual_frames.append(pd.read_csv(input_path))
all_data = pd.concat(individual_frames)

system_map = {}
for name, group in all_data.groupby("Cached(system_name)"):
    system_map[name] = group
    
for key in system_map:
    output_path = "%s/%s.txt" % (output_dir,key)
    make_bcl_file(output_path,system_map[key])

            
