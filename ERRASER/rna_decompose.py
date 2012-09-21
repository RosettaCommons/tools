#!/usr/bin/env phenix.python
import os.path
import imp

file_path = os.path.split( os.path.abspath(__file__) ) [0]
imp.load_source('erraser_util', file_path + '/erraser_util.py')

from erraser_util import *

input_pdb = sys.argv[1]

total_res = get_total_res(input_pdb)
n_chunk = int(total_res / 100.0 + 0.5)

if n_chunk <= 1 :
    print "Input pdb < 150 residues, no slicing is required."
else :
    print "Input pdb >= 150 residus, slice into %s chunks and minimize each one sequentially." % n_chunk
    res_slice_list = pdb_slice_into_chunks(input_pdb, n_chunk)
    current_chunk = 1
    for current_chunk, res_slice in enumerate(res_slice_list) :
        print "Chunk %s residues: %s" % (current_chunk+1, res_slice)

