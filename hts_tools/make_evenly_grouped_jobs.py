#!/usr/bin/env python 
import json
from optparse import OptionParser


def init_options():
    usage = "%prog --n_chunks input_screening_file.js output_prefix"
    parser=OptionParser(usage)
    parser.add_option("--n_chunks",dest="n_chunks",default=4)
    return parser
    
def get_emptiest_bin(bin_map):
    smallest_bin = 0
    smallest_bin_size = 100000
    for bin_id in bin_map:
        bin_size = 0
        for record in bin_map[bin_id]:
            bin_size += record[0]
        if smallest_bin_size > bin_size:
            smallest_bin = bin_id
            smallest_bin_size = bin_size
    return smallest_bin

if __name__ == "__main__":
    
    (options,args) = init_options().parse_args()
    job_file = args[0]
    prefix = args[1]
    bin_count = int(options.n_chunks)
    
    with open(job_file) as infile:
        jobs = json.load(infile)

    job_list = []
    all_jobs = 0
    for job in jobs:
        n_proteins = len(job["proteins"])
        n_ligands = len(job["ligands"])
        total_size = n_proteins*n_ligands
        job_list.append( (total_size,job) )
        all_jobs += total_size
    print all_jobs

    job_list = sorted(job_list,key=lambda x: x[0],reverse=True)

    bin_data = {}
    for bin_id in xrange(1,bin_count+1):
        bin_data[bin_id] = []

    for job_size,job_data in job_list:
        #print job_size, job_data
        bin_to_use = get_emptiest_bin(bin_data)
        bin_data[bin_to_use].append((job_size,job_data))

    for bin_id in bin_data:
        bin_size = 0
        for record in bin_data[bin_id]:
            bin_size += record[0]
        print bin_size

    for bin_id in bin_data:
        new_jobs = [record[1] for record in bin_data[bin_id]]
        with open("%s_%02d.js" % (prefix,bin_id),'w') as outfile:
            json.dump(new_jobs,outfile,indent=1)
