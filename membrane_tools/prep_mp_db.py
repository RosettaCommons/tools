#!/usr/bin/python
## @file prep_mpdb.py
##
## @brief 	Prep Data files for individual membrane proteins
## @details For a given membrane protein run, grabs tr pdb, splits by chain, 
##			generates a span file from structure, generates a lipid accessibility file, 
##			outputs a fasta file and generates and xml file based on the setup of these resources for
##			Rosetta membrane
##
## @author Rebecca Alford
## @note Last Modified: 1/26/14

# Python Headers
import sys, os
import argparse
import numpy
from subprocess import call

# Membrane Specific Import Headers
import get_transformed_pdb
import span_from_pdb
import lipid_acc_from_span
import write_mp_xml

## Write Embedding template (to be filled in by the user)
def write_embedding_template( chains ): 

	for chain_id in chains: 

		# Open writeable embed file
		filepath = chain_id + ".embed"
		f = open( filepath, 'w' )

		# Write template elements
		f.write( "POSE " + PDBID )
		f.write( "PROPERTY <XXX>" )
		f.write( "" )
		f.write( "CENTER X Y Z <METHOD_TAG>" )
		f.write( "NORMAL X Y Z <METHOD_TAG>" )
		f.write( "DEPTH D" )
		f.write( "" )
		f.write( "RESIDUES <XXX>" )


## Write Flags File for running integration tests (can be easily changed/commented out)
def write_flags(): 

	flags = "-run:constant_seed # Basic: Set a constant seed\n" \
 		"-run:jran 11111111 # Set constant seed to same as UTs\n"
 		"-overwrite # Basic: Allow overwrites for rerunning of unit tests\n" \
 		"-membrane:no_interpolate_Mpair # Membrane Scoring: Turn off interpolation betwene membrane layers \n" \
 		"-membrane:Menv_penalties # Membrane Scoring: Turn on membrane penalties for each run \n" \
		"-Membed_init # Membrane Scoring: Use initial embeddings (for now) \n" \
		"-in:file:membrane_chains chains.txt # Specify chains for PDB of interest\n"

	f = open( 'flags', 'w' ) 
	f.write(flags)

## Write Chain File
def write_chains( chains ): 

	# Print all Chains to a chains file
	g.open( "chains.txt", 'w' )
	for chain_id in chains: 
		g.write( chain_id + "\n")

## Write Rosetta command for running integration tests (can be changed/ commented out)
def write_cmd( PDBID ):

	cmd = "$rosetta_root/main/source/bin/mpframework_integration.macosclangdebug -database $rosetta_root/database "\
		  "@flags -jd2:resource_definition_files " + PDBID + ".xml"


## Given a PDB_ID, make a new directory and use this as the working directory
def prep_db( PDBID, symmetric, biomolecule ):

	# Make a new directory to store output files and then move to that dir
	os.mkdir( PDBID )
	os.chdir( PDBID )

	# Grab the transformed PDB
	get_transformed_pdb( PDBID )

	# Clean PDB using bcl and output corresponding fasta files
	if ( not symmetric ):
		os.system("$bcl PDBConvert " + pdb + "-bcl_pdb Split -output_prefix " + PDBID " -fasta ")
	else: 
		os.system("$bcl PDBConvert " + pdb + "-bcl_pdb Split -output_prefix " + PDBID " -fasta -biomolecule " + biomolecule )
	print pdb

	# Make Chain List from fasta sequences
	chains = []
	fasta_seqs = [ f for f in os.listdir("*.fasta") if os.path.isfile(f) ]
	for f in fasta_seqs: 
	     if ".fasta" in f:
	     	chains.append( os.path.splitext(f)[0] )



	## For each chain, get the appropriate span file and lipid acc data
	for chain_id in chains
		write_spanfile(chain_id)
		write_lipofile(chain_id)
		write_embedding_template(chain_id)

	## Create an xml file
	outxml = PDBID + ".xml"
	write_mp_xml( os.getcwd(), "chains.txt", outxml )

	## Write a flags file
	write_flags()

	## Create a rosetta commandfile to easily run from a cluster
	write_cmd()

# Print Help
def print_help():
    
    print "prep_mp_db.py --pdb <pdb> --mp_chains <chain ids> --non_mp_chains <chain ids> -symmetric <true/false> -biomolecule <##>"

    print "Required Args: "
    print "pdb = ID of the pdb of interest"
    print "mp_chains = Chains in the PDB that are considered membrane spanning"
    print "non_mp_chains = Chains in the PDB considered non membrane spanning"
    print "chain id = The chain id you are interested in. If more than one chain, "
    print "you cna pass the chain id without spaces. For example \"AB\" gets you, "
    print "chain A and B. \"A\" gets you chain A."
    print "\n",
    print "Optional Args: "
    print "biomolecue: Make the biomolecule <##> using the bcl"
    print "written by Rebecca Alford (Gray Lab)."
    sys.exit()

def main( argv ):

	# Parse Args
	parser = argparse.ArgumentParser( description='Generates required data inputs and resource definition XML for a given protein. Runs with the membrane framework')
	parser.add_argument('--pdb', help="PDB ID to prep")
	parser.add_argument('--mp_chains', help="Membrane chains")
	parser.add_argument('--non_mp_chains', help="Non membrane chains")
	parser.add_argument('-symmetric', help="Make a biomolecule given a single subunit in the PDB")
	parser.add_argument('-biomolecule', help="Biomolecule ## in the PDB")
	parser.parse_args()

	print "Writing input files for " + args.pdb

	# Prep PDB
	prep_db( args.pdb, args.symmetric, args.biomolecule )
	
	print "Done!"


if __name__ == "__main__" : main(sys.argv)
