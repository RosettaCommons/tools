#!/usr/bin/python

##########################################################

from numpy import *
import string
from make_tag import make_tag_with_dashes
from read_pdb import read_pdb

##########################################################

def get_surrounding_res( pdbfile, sample_res_list=[], radius=None, verbose=False ):

	### Make sure sample_res_list is a list of integers [ 4, 5, 6, 7 ]
	if ( len( sample_res_list ) == 1 ) and ( '-' in sample_res_list[0] ): 	
		sample_res_list = [ x for x in xrange( int(sample_res_list[0].split('-')[0]), int(sample_res_list[0].split('-')[-1])+1 ) ]
	else:	
		sample_res_list = [ int( res ) for res in sample_res_list ]
	
	### NEED TO TAKE CHAIN ID INTO ACCOUNT



	### Read PDB file
	( coords, pdb_lines, sequence, chains, residues ) = read_pdb( pdbfile )

	### Find surrounding residues if radius != None
	if radius:

		keep_res_list = []

		for surr_rsd in residues:
			is_surrounding_res = False		

			for sample_rsd in sample_res_list:
				if is_surrounding_res:	break

				for surr_atom in coords[chains[surr_rsd-1]][surr_rsd]:
					if is_surrounding_res:	break
					surr_atom_xyz = coords[chains[surr_rsd-1]][surr_rsd][surr_atom] 

					for sample_atom in coords[chains[sample_rsd-1]][sample_rsd]:
						if is_surrounding_res:	break
						sample_atom_xyz = coords[chains[sample_rsd-1]][sample_rsd][sample_atom]

						distance = sqrt( sum( power( subtract( surr_atom_xyz, sample_atom_xyz ) , 2 ) ) )
						if ( power( distance, 2 ) < power( float(radius), 2 ) ):
							if verbose:	print "res ", surr_rsd, " is a surrounding res, distance() = ", distance
							keep_res_list.append( surr_rsd )
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
	else:	keep_res_list = residues 

	### Remove sample residues from keep_res_list
	surrounding_res_list = [ res for res in keep_res_list if res not in sample_res_list ]
	return surrounding_res_list

##########################################################

def get_surrounding_res_tag( pdbfile, sample_res_list=[], radius=None, verbose=False ):

	surrounding_res_list = get_surrounding_res( pdbfile, sample_res_list=sample_res_list, radius=radius, verbose=verbose )

	if len( surrounding_res_list ):
		( coords, pdb_lines, sequence, chains, residues ) = read_pdb( pdbfile )
		surrounding_res_chains = [ chains[res-1] for res in surrounding_res_list ]
		surrounding_res_tag = make_tag_with_dashes( surrounding_res_list, surrounding_res_chains )
	else:	
		surrounding_res_tag = ''

	return surrounding_res_tag

##########################################################
##########################################################

if __name__=='__main__':

	import argparse

	parser = argparse.ArgumentParser(description='.')
	parser.add_argument('pdbfile')
	parser.add_argument('-sample_res', nargs='+', default=[])
	parser.add_argument('-radius', default=None)
	parser.add_argument('-make_tag', default=True)
	parser.add_argument('-verbose', default=False)
	args=parser.parse_args()

	if args.make_tag == True:
		surrounding_res_tag = get_surrounding_res_tag( args.pdbfile, sample_res_list=args.sample_res, radius=args.radius, verbose=args.verbose )
		print '\nSurrounding Residues: ',surrounding_res_tag.replace(' ',',')
	else: 
		surrounding_residues = get_surrounding_res( args.pdbfile, sample_res_list=args.sample_res, radius=args.radius, verbose=args.verbose )
		print '\nSurrounding Residues: ', string.join( [ str(x) for x in surrounding_residues ], ' ' ) 
