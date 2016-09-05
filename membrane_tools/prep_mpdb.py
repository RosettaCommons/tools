#!/usr/bin/env python
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

import sys, os
import argparse
import numpy
import shutil
from subprocess import call

import get_transformed_pdb
import span_from_pdb
import write_mp_xml

## @desc Prepare a pdb ID for integration testing and generate required mp 
##		 framework files
##
## @dependencies: clean_pdb, amino_acids, span_from_pdb, get_transformed_pdb, 
##				  write_mp_xml, remove_chain.pl
##
## Steps for Generating Data (used by the script)
##		(1) From PDB id, create new directory and make this your working directory
##		(2) Grab transformed PDB from the server (PDB TM Database)
##		(3) Generate one clean pdb per chain and a combined pbd
##		(4) Remove chains from the original PDB (this prevents making a copy of the PDB - more
##			memory efficient given that I currently do not need the whole PDB)
##		(5) Generate 1 spanfile per PDB chain 
##		(6) Generate a unified spanfile for the nonspecific PDB (maybe do this with a --hack)
##		(7) Generate an embedding file template for each PDB chain
##		(8) Write flags file
##	    (9) Write an XML resource definition file
##	   (10) Write a command file (and make it executable)
##	   (11) Done!

#===================================
## Environment variables (CHANGE THESE)
## Path to Rosetta executable
rosetta_exec = "/Users/rfalford12/apps/Rosetta/main/source/bin" # path to rosetta executables
rosetta_db = "/Users/rfalford12/apps/Rosetta/main/database" # path to rosetta database
remove_chain = "~/mpdocking/code/pdb_prep/remove_chain.pl" # path to remove chain.pl script
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

	# parse arguments
	args = parser.parse_args()

	print "Writing input files for ", args.pdb

	# Prep PDB
	prep_db( args.pdb, args.chains, "true" )

	print "Warning - you think you are done - but you still need to fill in data in the .embed files!"
	print "Done!"

#===================================
# Print Help
def print_help():
    
    print "prep_mp_db.py --pdb <pdb> --chains <chains to include>"
    print "Required Args: "
    print "pdb = ID of the pdb of interest"
    print "chains = chains to include in analysis"
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

	# Create a Cleaned pdb for each chain
	pdbpath = PDBID + "_tr"
	for i in range( len(chains) ): 
		clean_pdb_cmd = "python ../clean_pdb.py " + pdbpath + " " + chains[i]
		os.system( clean_pdb_cmd )

	# Create a Cleaned PDB for the full structure and remove chain IDs
	chainstring = ""
	for i in range( len(chains) ):
		chainstring = chainstring + chains[i]

	clean_pdb_cmd = "python ../clean_pdb.py " + pdbpath + " " + chainstring
	os.system( clean_pdb_cmd )
	
	result = PDBID.upper() + "_TR_" + chainstring + ".pdb"
	remove_chain_command = remove_chain + " -pdbfile " + result + " -outfile nochain.pdb"
	os.system( remove_chain_command )

	# Write chains.txt to dir
	write_chain_ids( PDBID, chains )

	# Build spanfile from whole pdb file
	span_from_pdb.get_spans_from_pdb( "nochain.pdb", "_", PDBID )

	## Create one span file per chain
	for c in chains: 
		
		my_id = PDBID.upper()
		my_chain = c.upper() 

		# Locate corresponding pdbfile
		pdbfile = my_id + "_TR_" + my_chain + ".pdb"
		newpath = my_id + "_" + my_chain + ".pdb"
		cpcmd = "cp " + pdbfile + " " newpath
		os.system( cpcmd )
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
