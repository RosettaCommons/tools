#!/usr/bin/env python

##########################################################

from numpy import *
import string
from make_tag import make_tag_with_dashes
from parse_tag import parse_tag
from read_pdb import read_pdb

##########################################################

def get_surrounding_res( pdbfile, sample_res_list=[], radius=None, verbose=False, keep_res=[] ):
	
	### Parse sample_res_list as a string
	sample_res_list, sample_chain_list = parse_tag( sample_res_list ) 
	
	### Read PDB file
	( coords, pdb_lines, sequence, chains, residues ) = read_pdb( pdbfile )

	### Find surrounding residues if radius != None
	if radius:


		if not keep_res is None:
			keep_res_list, keep_res_chain_list = parse_tag( keep_res )
		else:
			keep_res_list, keep_res_chain_list = [], []

		for surr_rsd_idx, surr_rsd in enumerate(residues):
			is_surrounding_res = False		

			for sample_rsd_idx, sample_rsd in enumerate(sample_res_list):
				if is_surrounding_res:	break

				for surr_atom in coords[chains[surr_rsd_idx]][surr_rsd]:
					if is_surrounding_res:	break
					surr_atom_xyz = coords[chains[surr_rsd_idx]][surr_rsd][surr_atom] 

					for sample_atom in coords[sample_chain_list[sample_rsd_idx]][sample_rsd]:
						if is_surrounding_res:	break
						sample_atom_xyz = coords[sample_chain_list[sample_rsd_idx]][sample_rsd][sample_atom]

						distance = sqrt( sum( power( subtract( surr_atom_xyz, sample_atom_xyz ) , 2 ) ) )
						if ( power( distance, 2 ) < power( float(radius), 2 ) ):
							if verbose:	print("res "+str(chains[surr_rsd_idx])+':'+str(surr_rsd)+" is a surrounding res, distance() = ", distance)
							keep_res_list.append( surr_rsd )
							keep_res_chain_list.append( chains[surr_rsd_idx] )
							is_surrounding_res = True
							break

		if not len( keep_res_list ):

			if not len( sample_res_list ):
				if verbose:	print('WARNING: len(sample_res_list) == ', len(sample_res_list), ' but radius == ', radius)
				if verbose:	print('Must supply sample_res to find residues within the given radius.' )
			else:
				if verbose:	print('WARNING: len(keep_res_list) == ', len(surrounding_res_list), ' for radius == ', radius)
				if verbose:	print('Try an expanded radius.')
		
	### All residues are surrounding if radius == None
	else:	
		keep_res_list = residues 
		keep_res_chain_list = chains

	assert( len(keep_res_list) == len(keep_res_chain_list) )

	### Remove sample residues from keep_res_list
	surrounding_res_list = []
	surrounding_res_chain_list = []
	for rsd_idx, rsd in enumerate(keep_res_list):
		chain = keep_res_chain_list[ rsd_idx ]
		if rsd in sample_res_list:
			if chain == sample_chain_list[ sample_res_list.index(rsd) ]:  continue
		surrounding_res_list.append( rsd )
		surrounding_res_chain_list.append( chain )

	assert( len(surrounding_res_list) == len(surrounding_res_chain_list) )

	return surrounding_res_list, surrounding_res_chain_list

##########################################################

def get_ellipsoid_envelope_res( pdbfile, sample_res_list=[], radius=None, verbose=False ):
	

	### Parse sample_res_list as a string
	sample_res_list, sample_chain_list = parse_tag( sample_res_list ) 

	### Read PDB file
	( coords, pdb_lines, sequence, chains, residues ) = read_pdb( pdbfile )

	'''
	five_prime_boundary  = sample_res_list[0] - 1
	three_prime_boundary = sample_res_list[-1] + 1
	num_loop_res	 = three_prime_boundary - five_prime_boundary

	if verbose:
		print('five_prime_boundary  = '+str(five_prime_boundary))
		print('three_prime_boundary = '+str(three_prime_boundary))
		print('num_loop_res		    = '+str(num_loop_res))

	############ PORT TO PYTHON
	five_prime_foci  = pose.residue( five_prime_boundary ).xyz( " O3'" )
	three_prime_foci = pose.residue( three_prime_boundary ).xyz( " C5'" )
	############


	expand_radius = float(radius)
	int_expand_radius = int( 10.*expand_radius)

	foci_sep_dist = sqrt( sum( power( subtract( three_prime_foci, five_prime_foci ), 2 ) ) )

	O3I_C5I_PLUS_ONE_MAX_DIST = 3.968000
	O3I_O3I_PLUS_ONE_MAX_DIST = 7.45583
	max_loop_length = ( num_loop_res*O3I_O3I_PLUS_ONE_MAX_DIST ) + O3I_C5I_PLUS_ONE_MAX_DIST

	major_diameter = max_loop_length + expand_radius

	if verbose:
		print("foci_sep_dist   = "+str(foci_sep_dist))
		print("major_diameter  = "+str(major_diameter))
		print("max_loop_length = "+str(max_loop_length))
		print("expand_radius   = "+str(expand_radius))




	############ Create a rotated pose which have the major axis as the z-axis and the center of the ellipsoid as the origin!/////////
	new_origin = ( five_prime_foci + three_prime_foci ) / 2.0
	
	internal_z_axis = [ 0.0, 0.0, 1.0 ]

	z_axis = ( three_prime_foci - five_prime_foci )
	z_axis = np.normalize( z_axis )

	y_axis = cross( internal_z_axis, z_axis )
	y_axis = np.normalize( y_axis )

	x_axis = cross( y_axis, z_axis )
	x_axis = np.normalize( x_axis )

	coordinate_matrix = [ x_axis, y_axis, z_axis ];
	rotation_matrix = np.inverse( coordinate_matrix );

		for ( Size i = 1; i <= rotated_pose.total_residue(); ++i ) {
		for ( Size j = 1; j <= rotated_pose.residue_type( i ).natoms(); ++j ) { // use residue_type to prevent internal coord update

			id::AtomID const id( j, i );

			rotated_pose.set_xyz( id, rotated_pose.xyz( id ) - new_origin );
			rotated_pose.set_xyz( id, rotation_matrix * rotated_pose.xyz( id ) );

		}
	}

	dump_pdb( rotated_pose, "rotated_" + pose_name );

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ){

		if ( sample_res_list.has_value( seq_num ) ){
			std::cout << "res " << seq_num << " is a sample_res" << std::endl;
			keep_res_list.push_back( seq_num );
			continue;
		}

		if ( additional_slice_res_list.has_value( seq_num ) ){
			std::cout << "res " << seq_num << " is in additional_slice_res_list" << std::endl;
			keep_res_list.push_back( seq_num );
			continue;
		}


		core::conformation::Residue const & surrounding_rsd = pose.residue( seq_num );

		core::conformation::Residue const & rotated_surr_rsd = rotated_pose.residue( seq_num );

		if ( surrounding_rsd.natoms() !=	rotated_surr_rsd.natoms()	){
			utility_exit_with_message( "surrounding_rsd.natoms() !=	rotated_surr_rsd.natoms()" );
		}

		for ( Size atomno = 1; atomno <= surrounding_rsd.natoms(); atomno++ ){

			/////////////////////////////////////////////////////////////////////////////

			numeric::xyzVector< Real > const atom_xyz = surrounding_rsd.xyz( atomno );

			Real const sum_distances_to_foci = ( atom_xyz - five_prime_foci ).length() + ( atom_xyz - three_prime_foci ).length();

			bool const method_1 = ( sum_distances_to_foci <= major_diameter ) ? true : false;

			//////////////////////////////////////////////////////////////////////////////
			numeric::xyzVector< Real > const rotated_atom_xyz = rotated_surr_rsd.xyz( atomno );

			Real const r_sq = rotated_atom_xyz.x()*rotated_atom_xyz.x() + rotated_atom_xyz.y()*rotated_atom_xyz.y();
			Real const z_sq = rotated_atom_xyz.z()*rotated_atom_xyz.z();

			Real const denominator_one = ( ( major_diameter*major_diameter ) - ( foci_sep_dist*foci_sep_dist ) )/4.0;

			Real const denominator_two = ( major_diameter*major_diameter )/4.0;

			Real const ellipse_equation = ( r_sq/denominator_one ) + ( z_sq/denominator_two );

			bool const method_2 = ( ellipse_equation <= 1.0 ) ? true: false;

			//////////////////////////////////////////////////////////////////////////////


			if ( method_1 != method_2 ){
				std::cout << "method_1 != method_2" << std::endl;
				output_boolean( "method_1 = ", method_1, TR ); std::cout << std::endl;
				output_boolean( "method_2 = ", method_2, TR ); std::cout << std::endl;
				utility_exit_with_message( "method_1 != method_2" );
			}


			if ( method_1 ){
				if ( keep_res_list.has_value( seq_num ) == false ){ //Not already in the list!
					std::cout << "seq_num ( " << seq_num << " ) is a surrounding_res, sum_distances_to_foci = " << sum_distances_to_foci << " major_diameter = " << major_diameter << std::endl;
					keep_res_list.push_back( seq_num );
					//break;
				}
			}
		}
	}

	'''

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
							if verbose:	print("res "+str(chains[surr_rsd_idx])+':'+str(surr_rsd)+" is a surrounding res, distance() = ", distance)
							keep_res_list.append( surr_rsd )
							keep_res_chain_list.append( chains[surr_rsd_idx] )
							is_surrounding_res = True
							break

		if not len( keep_res_list ):

			if not len( sample_res_list ):
				if verbose:	print('WARNING: len(sample_res_list) == ', len(sample_res_list), ' but radius == ', radius)
				if verbose:	print('Must supply sample_res to find residues within the given radius.' )
			else:
				if verbose:	print('WARNING: len(keep_res_list) == ', len(surrounding_res_list), ' for radius == ', radius)
				if verbose:	print('Try an expanded radius.')
		
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

def get_surrounding_res_tag( pdbfile, sample_res_list=[], radius=None, verbose=False, csv=False, ellipsoid_envelope_mode=False, keep_res=None ):

	if not ellipsoid_envelope_mode:
		( surrounding_res_list, surrounding_res_chain_list ) = get_surrounding_res( pdbfile, sample_res_list=sample_res_list, radius=radius, verbose=verbose, keep_res=keep_res )
	else:
		( surrounding_res_list, surrounding_res_chain_list ) = get_ellipsoid_envelope_res( pdbfile, sample_res_list=sample_res_list, radius=radius, verbose=verbose )

	if len( surrounding_res_list ):		surrounding_res_tag = make_tag_with_dashes( surrounding_res_list, surrounding_res_chain_list )
	else:								surrounding_res_tag = ''
	if csv:								surrounding_res_tag = surrounding_res_tag.replace( ' ', ',' )
	if surrounding_res_tag[0] == ',':	surrounding_res_tag = surrounding_res_tag[1:]
	return surrounding_res_tag

##########################################################
##########################################################

if __name__=='__main__':

	import argparse

	parser = argparse.ArgumentParser(description='.')
	parser.add_argument('pdbfile')
	parser.add_argument('-sr','--sample_res', nargs='+', default=[])
	parser.add_argument('-r','--radius', default=None)
	parser.add_argument('--ellipsoid_envelope_mode', action="store_true")
	parser.add_argument('--make_tag', action="store_true")
	parser.add_argument('--make_tag_csv', action="store_true")
	parser.add_argument('-v', '--verbose', help="increase output verbosity", action="store_true")	
	args=parser.parse_args()


	surrounding_res_tag = get_surrounding_res_tag( args.pdbfile, sample_res_list=args.sample_res, radius=args.radius, verbose=args.verbose, csv=args.make_tag_csv, ellipsoid_envelope_mode=args.ellipsoid_envelope_mode )

	if not args.make_tag and not args.make_tag_csv:	
		surrounding_residues, surrounding_chains = parse_tag( surrounding_res_tag )
		surrounding_res_tag = string.join( [ str(x) for x in surrounding_residues ], ' ' )

	if args.verbose:	print('\nSurrounding Residues: ', surrounding_res_tag)#.replace(' ',',')
	else:				print(surrounding_res_tag)
	
