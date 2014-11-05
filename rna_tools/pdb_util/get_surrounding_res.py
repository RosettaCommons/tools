#!/usr/bin/python

##########################################################

from numpy import *
import string
from make_tag import make_tag_with_dashes
from parse_tag import parse_tag
from read_pdb import read_pdb

##########################################################

def get_surrounding_res( pdbfile, sample_res_list=[], radius=None, verbose=False ):
	
	### Parse sample_res_list as a string
	sample_res_list, sample_chain_list = parse_tag( sample_res_list ) 

	### Read PDB file
	( coords, pdb_lines, sequence, chains, residues ) = read_pdb( pdbfile )

	### Find surrounding residues if radius != None
	if radius:

		keep_res_list = []
		keep_res_chain_list = []


		for surr_rsd_idx in xrange( len(residues) ):
			is_surrounding_res = False		
			surr_rsd = residues[ surr_rsd_idx ]

			for sample_rsd_idx in xrange( len(sample_res_list) ):
				if is_surrounding_res:	break
				sample_rsd = sample_res_list[ sample_rsd_idx ]

				for surr_atom in coords[chains[surr_rsd_idx]][surr_rsd]:
					if is_surrounding_res:	break
					surr_atom_xyz = coords[chains[surr_rsd_idx]][surr_rsd][surr_atom] 

					for sample_atom in coords[sample_chain_list[sample_rsd_idx]][sample_rsd]:
						if is_surrounding_res:	break
						sample_atom_xyz = coords[sample_chain_list[sample_rsd_idx]][sample_rsd][sample_atom]

						distance = sqrt( sum( power( subtract( surr_atom_xyz, sample_atom_xyz ) , 2 ) ) )
						if ( power( distance, 2 ) < power( float(radius), 2 ) ):
							if verbose:	print "res ", chains[surr_rsd_idx],':',surr_rsd, " is a surrounding res, distance() = ", distance
							keep_res_list.append( surr_rsd )
							keep_res_chain_list.append( chains[surr_rsd_idx] )
							is_surrounding_res = True
							break

		if not len( keep_res_list ):

			if not len( sample_res_list ):
				if verbose:	print 'WARNING: len(sample_res_list) == ', len(sample_res_list), ' but radius == ', radius
				if verbose:	print 'Must supply sample_res to find residues within the given radius.' 
			else:
				if verbose:	print 'WARNING: len(keep_res_list) == ', len(surrounding_res_list), ' for radius == ', radius
				if verbose:	print 'Try an expanded radius.'
		
	### All residues are surrounding if radius == None
	else:	
		keep_res_list = residues 
		keep_res_chain_list = chains

	assert( len(keep_res_list) == len(keep_res_chain_list) )

	### Remove sample residues from keep_res_list
	surrounding_res_list = []
	surrounding_res_chain_list = []
	for rsd_idx in xrange( len( keep_res_list ) ):
		rsd = keep_res_list[ rsd_idx ]
		chain = keep_res_chain_list[ rsd_idx ]
		if rsd in sample_res_list: continue
		surrounding_res_list.append( rsd )
		surrounding_res_chain_list.append( chain )

	assert( len(surrounding_res_list) == len(surrounding_res_chain_list) )

	return surrounding_res_list, surrounding_res_chain_list

##########################################################

def get_surrounding_res_tag( pdbfile, sample_res_list=[], radius=None, verbose=False ):
	surrounding_res_list, surrounding_res_chain_list = get_surrounding_res( pdbfile, sample_res_list=sample_res_list, radius=radius, verbose=verbose )
	if len( surrounding_res_list ):	surrounding_res_tag = make_tag_with_dashes( surrounding_res_list, surrounding_res_chain_list )
	else:							surrounding_res_tag = ''
	return surrounding_res_tag

##########################################################
##########################################################

if __name__=='__main__':

	import argparse

	parser = argparse.ArgumentParser(description='.')
	parser.add_argument('pdbfile')
	parser.add_argument('-sr','--sample_res', nargs='+', default=[])
	parser.add_argument('-r','--radius', default=None)
	parser.add_argument('--make_tag', action="store_true")
	parser.add_argument('--make_tag_csv', action="store_true")
	parser.add_argument('-v', '--verbose', help="increase output verbosity", action="store_true")	
	args=parser.parse_args()


	surrounding_res_tag = get_surrounding_res_tag( args.pdbfile, sample_res_list=args.sample_res, radius=args.radius, verbose=args.verbose )
	surrounding_residues, surrounding_chains = parse_tag( surrounding_res_tag )
		
	if args.make_tag_csv:	
		surrounding_res_tag = surrounding_res_tag.replace( ' ', ',' )
		if surrounding_res_tag[0] == ',': surrounding_res_tag = surrounding_res_tag[1:]
	elif not args.make_tag:	surrounding_res_tag = string.join( [ str(x) for x in surrounding_residues ], ' ' )

	if args.verbose:	print '\nSurrounding Residues: ', surrounding_res_tag#.replace(' ',',')
	else:				print surrounding_res_tag
	