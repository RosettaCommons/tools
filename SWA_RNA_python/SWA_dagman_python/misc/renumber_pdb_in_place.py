#!/usr/bin/env python

import string
from sys import argv,stderr,stdout
from os import popen,system
from os.path import exists

from SWA_dagman_python.database.SWA_amino_acids import longer_names
from SWA_dagman_python.utility.PDB_operations import renumber_pdb_in_place_func


assert( len(argv)==2)

pdbnames = argv[1:]

if(len(pdbnames)!=1): error_exit_with_message("len(pdbnames)!=1")

renumber_pdb_in_place_func(pdbnames[0])
