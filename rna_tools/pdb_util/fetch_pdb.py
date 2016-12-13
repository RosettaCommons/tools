#!/usr/bin/env python

import urllib
from os import getcwd
from os.path import exists
import argparse

#######################################################################################

def fetch(pdb_id):
	url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdb_id.upper()
	return urllib.urlopen( url ).read()

def fetch_pdb(pdb_id, inpath_dir=getcwd()):
	if inpath_dir.split()[-1] != '/': inpath_dir+= '/'
	pdb_out = inpath_dir+pdb_id+'.pdb'
	if not exists( pdb_out ):	open( pdb_out, 'w' ).write( fetch( pdb_id ) )
	return pdb_out

#######################################################################################

if __name__=='__main__':
	
	parser = argparse.ArgumentParser(description='Download PDB files from RCSB.')
	parser.add_argument('pdb_id', help='PDB ID to fetch')
	parser.add_argument('-inpath_dir', help='path to the directory in which the PDB should be written', default=getcwd())
	args = parser.parse_args()

	fetch_pdb(args.pdb_id, inpath_dir=args.inpath_dir)