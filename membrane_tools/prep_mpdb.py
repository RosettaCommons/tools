#!/usr/bin/python
## @file prep_mpdb.py
##
## @brief 	Prep Data files for individual membrane proteins
## @details For a given membrane protein run, grabs tr pdb, splits by chain, 
##			generates a span file from structure, generates a lipid accessibility file, 
##			outputs a fasta file and generates and xml file based on the setup of these 
##			resources for Rosetta membrane
##
## @author Rebecca Alford
## @note Last Modified: 2/5/14

# Python Headers
import sys, os
import argparse
import numpy
import shutil
from subprocess import call

# Membrane Specific Import Headers
import get_transformed_pdb
import span_from_pdb
import write_mp_xml

# Notes from 2/6/14 - removing bcl dependency. trpdbs should already have 
# symmetric molecules in them. Clean_pdb is a much better way to do the other thing...

## @desc Prepare a pdb ID for integration testing and generate required mp 
##		 framework files
##
## @dependencies: clean_pdb, amino_acids, span_from_pdb, get_transformed_pdb, 
##				  write_mp_xml, bcl (for chain splitting, until I find/write an alternative)
##
## Steps for Generating Data: 
##		(1) From PDB id, create new directory and make this your working directory
##		(2) Grab transformed PDB from the server (PDB TM Database)
##		(3) From the TR_PDB, generate a clean pdb with the clean_pdb.py script and 
##			and chains desired from the user (fyi - will split the fasta but no pdb)
##		(4) Split the cleaned PDB by chain (no renumbering - but convert to natural aa type)
##		(5) Remove chains from the original PDB (this prevents making a copy of the PDB - more
##			memory efficient given that I currently do not need the whole PDB)
##		(6) Generate 1 spanfile per PDB chain 
##		(7) Generate a unified spanfile for the nonspecific PDB (maybe do this with a --hack)
##		(8) Generate an embedding file template for each PDB chain
##		(9) Write flags file
##	   (10) Write an XML resource definition file
##	   (11) Write a command file (and make it executable)
##	   (12) Done!

#===================================
## Environment variables
## Path to Rosetta executable
rosetta_exec = "/Users/rfalford12/apps/Rosetta/main/source/bin" #insert your path here
rosetta_db = "/Users/rfalford12/apps/Rosetta/main/database"
bcl = "/Users/rfalford12/apps/BioChemicalLibrary-3.0.0_mac/bcl" #insert your path here
## end of env vars

#===================================
# Main
def main( argv ):
	''' Script for preparing a PDB for integration testing by cleaning required PDB files
		and creating the necessary inputs to the Membrane Framework (RosettaMembrane)
	'''

	# Set up arg parsing
	parser = argparse.ArgumentParser( description='collecting all required components for running things with the membrane framework (RosettaMembrane)')

	# Load PDB File
	parser.add_argument('--pdb', '-f',
    action="store",
    help="Name of pdb filename. Should only contain a single model! Example: 'file.pdb'.", )

	# Specify Desired Chains
	parser.add_argument('--chains', '-c',
    action="store",
    nargs='+',
    help="Chain to extract from the pdb. Should be uppercase [A-Z].", )

    # Hack the clean pdb process to build full non-chain dependent pdb
   # parser.add_argument('--hack', '-k',
   # action="store",
   # help="Hack clean pdb to make an addiitonal full spanfile for scoring.", )

	#parse arguments
	args = parser.parse_args()

	print "Writing input files for ", args.pdb

	# Prep PDB
	prep_db( args.pdb, args.chains, "true" )
	
	print "Done!"

#===================================
# Print Help
def print_help():
    
    print "prep_mp_db.py --pdb <pdb> --chains <chains to include> --hack <true/false>"
    print "Required Args: "
    print "pdb = ID of the pdb of interest"
    print "chains = chains to include in analysis"
    print "hack: build full spanfile with renumbering for scoring in pre-iter version 2014"
    print "written by Rebecca Alford (Gray Lab)."
    sys.exit()

#===================================
## Given a PDB_ID, make a new directory and use this as the working directory
def prep_db( PDBID, chains, hack ):

	# Make a new directory to store output files and then move to that dir
	os.mkdir( PDBID )
	os.chdir( PDBID )

	# Grab the transformed PDB
	get_transformed_pdb.get_transformed_pdb( PDBID )

	# Check file exists
	if ( not os.path.isfile( PDBID + "_tr.pdb") ):
		sys.exit("Failed to download tr pdb from the server!")

	# Clean the PDB using clean_pdb script (os.sys becasue clean_pdb is not an object)
	pdbpath = PDBID + "_tr"
	chainstring = ""
	for i in range( len(chains) ):
		chainstring = chainstring + chains[i]

	clean_pdb_cmd = "python ../clean_pdb.py " + pdbpath + " " + chainstring
	print "Cleaning the PDB"
	os.system( clean_pdb_cmd )

	# hack - right now - becasue we have no lipid file dependency, we won't need the fastas
	# generated from clean pdb. However, in the future if you get that to work, it might be
	# a better idea to delete these and regenerate them with the bcl which is slightly
	# more sophisticated

	# Split existing PDB into multiple chains with additional cleaning steps
	clean_pdb_result = PDBID.upper() + "_TR_" + chainstring + ".pdb"
	bcl_command = bcl + " protein:PDBConvert " + clean_pdb_result + " -bcl_pdb Split -output_prefix " + PDBID + "_"
	os.system( bcl_command )

	# Write chains.txt to dir
	write_chain_ids( PDBID, chains )

	## Write non split/hacky span file
	if ( hack == "true" ):

		# Create a copy of the pdb with removed chains
		nochain_path = "nochain.pdb" ## change this...
		remove_chain_command = "~/mpdocking/code/pdb_prep/remove_chain.pl -pdbfile " + clean_pdb_result + " -outfile " + nochain_path
		os.system( remove_chain_command )

		# BUild spanfile from that pdb
		span_from_pdb.get_spans_from_pdb( nochain_path, "_", PDBID )

	## Create one span file per chain
	for c in chains: 
		
		my_id = PDBID.upper()
		my_chain = c.upper() 

		# Locate corresponding pdbfile
		pdbfile = my_id + "_" + my_chain + ".pdb"
		prefix = pdbfile[0:4] + "_"
		prefix2 = prefix + my_chain

		# Create Span File and write embedding template
		span_from_pdb.get_spans_from_pdb( pdbfile, my_chain, prefix )
		write_embedding_template( prefix2 )

	## Write flags file
	write_flags()

	## Write xml file
	outxml = PDBID + ".xml"
	write_mp_xml.setup_xml( os.getcwd(), "chains.txt", outxml )

	## Write final command path
	write_cmd( PDBID )

	## return
	return

########## File writing methods ########## 

#===================================
## Write chain IDs to a file
def write_chain_ids( PDBID, chains ): 

	# Make Chain Info List 
	# (much less hacky than the previous approach using fasta files)
	chain_list = []
	for c in chains: 
		my_id = PDBID.upper()
		my_chain = c.upper()
		chain_id = my_id + "_" + my_chain
		chain_list.append( chain_id )

	with open('chains.txt', 'a') as f:	
		for i in range( len(chain_list) ): 	
			f.write(chain_list[i] + '\n')

#===================================
## Write Embedding template (to be filled in by the user)
def write_embedding_template( chain_id ): 

	# Open writeable embed file
	filepath = chain_id + ".embed"
	with open( filepath, 'w' ) as f:

		# Write template elements
		f.write( "POSE " + chain_id + '\n' )
		f.write( "" + '\n' )
		f.write( "CENTER X Y Z <METHOD_TAG>" + '\n' )
		f.write( "NORMAL X Y Z <METHOD_TAG>" + '\n' )
		f.write( "DEPTH D" + '\n' )

		# be nice to the buffer
		f.close()

	# Check that the new file exists - otherwise throw an error
	if ( not os.path.isfile( filepath ) ):
		sys.exit("Embedding file template was incorrectly created" )

#===================================
## Write Rosetta command for running integration tests (can be changed/ commented out)
def write_cmd( PDBID ):

	# Info for command file
	cmdfile = PDBID + "_integration.sh"
	cmd = rosetta_exec + "/mpframework_integration.macosclangdebug -database " + rosetta_db + " @flags -jd2:resource_definition_files " + PDBID + ".xml"

	# Useful commands (at 35000ft, can't look up the os.sys equivalent for this)
	exec_cmd = "chmod +x " + cmdfile
	with open ( cmdfile, 'w') as f: 
		f.write( "#!/bin/bash" + '\n')
		f.write( cmd + '\n' )

	# Make result file an executable
	os.system( exec_cmd )

#===================================
## Write Flags File for running integration tests (can be easily changed/commented out)
def write_flags(): 

	# Define current flags set (can interchange this method if you want to define your own flags
	# or you don't need these - some are optional, some highly coupled to the integration test)
	flags = "-run:constant_seed # Basic: Set a constant seed\n" \
 		"-run:jran 11111111 # Set constant seed to same as UTs\n"\
 		"-overwrite # Basic: Allow overwrites for rerunning of unit tests\n" \
 		"-membrane:no_interpolate_Mpair # Membrane Scoring: Turn off interpolation betwene membrane layers \n" \
 		"-membrane:Menv_penalties # Membrane Scoring: Turn on membrane penalties for each run \n" \
		"-Membed_init # Membrane Scoring: Use initial embeddings (for now) \n" \
		"-in:file:membrane_chains chains.txt # Specify chains for PDB of interest\n"

	# Write data to file
	with open ( 'flags', 'w' ) as f: 
		f.write(flags)
		f.close()

	# Check that my file exists in the correct dir
	# losts of arbitrary files named flags would really suck....
	if ( not os.path.isfile( 'flags' ) ): 
		sys.exit("Missing required flags file!")

## end of file writing methods

if __name__ == "__main__" : main(sys.argv)
