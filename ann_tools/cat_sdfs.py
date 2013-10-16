#!/usr/bin/env python

import gzip
import glob
import sys
import os 

input_dir = sys.argv[1]
output_path = sys.argv[2]

outfile = gzip.open(output_path,'w')
for path in glob.glob(input_dir+"/*.sdf.gz"):
    infile = gzip.open(path)
    for line in infile:
        outfile.write(line)
    infile.close()
    
    os.unlink(path)
        