#!/usr/bin/env python3

import os
import argparse
import json
import pandas as pd 
import numpy as np
from array import array

# Import gremlin with modifications that return a raw tensor
import gremlin

# Dumps the raw GREMLIN tensor in the Rosetta MathNTensor compatible format.
# Authors: Moritz Ertelt (moritz.ertelt@stud.uni-regensburg.de), Samuel Schmitz (Samuel.Schmitz@vanderbilt.edu)

# Creating the ArgParser
pars = argparse.ArgumentParser()

pars.add_argument("MSA", help="The msa file to be analysed with GREMLIN")
pars.add_argument("OUTPUT", help="Name of the output files. Suffixes .bin, .index and .json will be appended")

args = pars.parse_args()

headers, seqs = gremlin.parse_fasta(str(args.MSA))

msa = gremlin.mk_msa(seqs, gap_cutoff=0.75, eff_cutoff=0.8)
mrf = gremlin.GREMLIN(msa, lam_w = 0.01)
mtx_single, mtx, apcTensor, apcTensor0 = gremlin.get_mtx(mrf)  

# Original Raw Tensor
# saving the Tensor and JSON file with shape information
output_file = open('%s.bin' % args.OUTPUT, 'wb')
float_arrayRaw = array('d', list(mrf["w"].flatten()))
float_arrayRaw.tofile(output_file)
output_file.close()

data = {}
data['n_bins'] = [mrf['w'].shape[0], 20, mrf['w'].shape[0], 20]
data['type'] = ['double']

with open('%s.json' % args.OUTPUT, 'w') as outfile:
    json.dump(data, outfile)

# Saving a index file to map the right positions to the tensor
np.savetxt("%s.index" % args.OUTPUT, mrf['v_idx'].astype(int),fmt='%i')
