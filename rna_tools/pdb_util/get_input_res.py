#!/usr/bin/python

##########################################################

import numpy as np
import string

import make_tag
from read_pdb import read_pdb

##########################################################

def distance( a, b ):
 	return np.sqrt( np.sum( np.power( np.subtract( a, b ) , 2 ) ) )

##########################################################

def get_peripheral_res( pdbfile, sample_res_list, radius ):

	( coords, pdb_lines, sequence, chains, residues ) = read_pdb( pdbfile )

	keep_res_list = []

	for seq_num in residues:
		is_surrounding_res = False		
		surr_rsd = seq_num

		for sample_res in sample_res_list:
			if is_surrounding_res:	break
			sample_rsd = sample_res

			for surr_atom in coords[chains[surr_rsd-1]][surr_rsd]:
				if is_surrounding_res:	break
				surr_atom_xyz = coords[chains[surr_rsd-1]][surr_rsd][surr_atom] 

				for sample_atom in coords[chains[sample_rsd-1]][sample_rsd]:
					if is_surrounding_res:	break
					sample_atom_xyz = coords[chains[sample_rsd-1]][sample_rsd][sample_atom]

					if ( np.power( distance( surr_atom_xyz, sample_atom_xyz ), 2 ) < np.power( radius, 2 ) ):
						print "res ", seq_num, " is a surrounding res, distance() = ", distance( surr_atom_xyz, sample_atom_xyz )
						keep_res_list.append( seq_num )
						is_surrounding_res = True
						break

	peripheral_res_list = [ res for res in residues if res not in keep_res_list and res not in sample_res_list ]

	return peripheral_res_list

##########################################################

def get_input_res_tag( pdbfile, sample_res_list=[], radius=None ):

	( coords, pdb_lines, sequence, chains, residues ) = read_pdb( pdbfile )

	input_res_list = get_input_res( pdbfile, sample_res_list=sample_res_list, radius=radius )
	
	if len(input_res_list): 
		input_chains = [ chains[res-1] for res in input_res_list ]
		input_res_tag = make_tag.make_tag_with_dashes( input_res_list, input_chains )
	else:
		print 'ERROR: len(input_res_list) == ', len(input_res_list)
		print 'Try a larger radius.'
		input_res_tag = ''

	return input_res_tag

##########################################################

def get_input_res( pdbfile, sample_res_list=[], radius=None ):
	
	if '-' in sample_res_list: sample_res_list = [ str(x) for x in xrange( int(sample_res_list.split('-')[0]), int(sample_res_list.split('-')[-1])+1 ) ]
	sample_res_list = [ int( res ) for res in sample_res_list ]

	( coords, pdb_lines, sequence, chains, residues ) = read_pdb( pdbfile )
	
	peripheral_res_list = []
	if radius:	peripheral_res_list = get_peripheral_res( pdbfile, sample_res_list, radius )
	
	input_res_list = [ res for res in residues if res not in sample_res_list and res not in peripheral_res_list ]

	return input_res_list

##########################################################
##########################################################

if __name__=='__main__':

	import argparse

	parser = argparse.ArgumentParser(description='.')
	parser.add_argument('pdbfile')
	parser.add_argument('-sample_res', nargs='+', default=[])
	parser.add_argument('-radius', default=None)
	args=parser.parse_args()

		
	input_res_tag = get_input_res_tag( pdbfile=args.pdbfile, sample_res_list=args.sample_res, radius=float(args.radius) )
	print 'Input Residues: ',input_res_tag
	