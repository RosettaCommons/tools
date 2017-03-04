#!/usr/bin/env python

from sys import stdout,argv
from os import system

from SWA_dagman_python.utility.PDB_operations import replace_chain_inplace_func


assert( len(argv)==3)

actualpdbnames = argv[1:-1]

if(len(actualpdbnames)!=1): error_exit_with_message("len(actualpdbnames)!=1")

newchain = argv[-1]

replace_chain_inplace_func(actualpdbnames[0], newchain)


