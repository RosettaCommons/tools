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
## @note Last Modified: 1/17/14

def main( argv ):


	# Parse Args
	parser = 
	# Parse Args
	parser = argparse.ArgumentParser( description='Generate a resource definition XML file given a list of membrane chains')
	parser.add_argument('--workdir', help="working directory for run - should be a full path from Rosetta/source")
	parser.add_argument('--chains', help="list of membrane chains to include in the run")
	parser.add_argument('--outfile', help='output xml file')
	parser.add_argument('--include_lips', help='include lipophobicity score in run')
	args = parser.parse_args()

	# Setup xml file
	setup_xml( args.workdir, args.chains, args.outfile, args.include_lips )

if __name__ == "__main__" : main(sys.argv)
