#!/usr/bin/python
## @file prep_membrane.py
##
## @brief Make membrane resource xml script
## @details Generate Membrane Resource Definition File from Given List of Chains
##
## @author Rebecca Alford
## @note Last Modified: 1/16/14

import sys, os
import subprocess
import argparse

# Create XML Template to be filled in by chain
def setup_xml( workdir, chains, outfile, lips_exp ):

	# Read membrane chains into an array
	mp_chains = []
	with open( chains, 'r') as mmembrane_chains:
		for line in mmembrane_chains:
			mp_chains.append(line.rstrip('\r\n'))

	# Open output xml file
	f = open( outfile, 'w' )

	# Write out xml header (will add some description headers later)
	f.write('<JD2ResourceManagerJobInputter>\n')
	f.write('    <ResourceLocators>\n')
	f.write('        <FileSystemResourceLocator tag="spanning_locator" base_path="' + workdir + '"/>\n')
	f.write('        <FileSystemResourceLocator tag="embedding_locator" base_path="' + workdir + '"/>\n')
	f.write('        <FileSystemResourceLocator tag="lipsexp_locator" base_path="' + workdir + '"/>\n')
	f.write('        <FileSystemResourceLocator tag="startstruct_locator" base_path="' + workdir + '"/>\n')
	f.write('    </ResourceLocators>')

	# Specifying a default set of resource options - can always add more!
	f.write('	    <ResourceOptions>\n')
	f.write('        <PoseFromPDBOptions\n')
	f.write('           tag="pdb1"\n')
	f.write('           ignore_unrecognized_res=1\n')
	f.write('        />\n')
	f.write('    </ResourceOptions>\n')

	# write out xml resource definitions for each chain
	f.write('    <Resources>\n')
	for i in range( len( mp_chains ) ):
		f.write('        <PoseFromPDB tag="' + mp_chains[i] + '_my" locator="startstruct_locator" locatorID="' + workdir + mp_chains[i] + '.pdb"/>\n')
		f.write('        <EmbedDef tag="' + mp_chains[i] + '_my_embed" locator="embedding_locator" locatorID="' + workdir + mp_chains[i] + '.embed"/>\n')
		f.write('        <SpanFile tag="' + mp_chains[i] + '_my_span" locator="spanning_locator" locatorID="' + workdir + mp_chains[i] + '.span"/>\n')
	f.write('    </Resources>\n')

	# write out xml job definition
	f.write('    <Jobs>\n')
	f.write('        <Job name="membrane">\n')
	f.write('            <Data desc="startstruct" resource_tag="' + mp_chains[i] + '_"/>\n') ### this is a pre-jd2 hack

	# Write resources per chain
	for j in range ( len( mp_chains ) ):
		f.write('			 <Data desc="' + mp_chains[i] + '_lips" resource_tag="' + mp_chains[i] + '_lips"/>\n')
		f.write('            <Data desc="' + mp_chains[i] + '_embed" resource_tag="' + mp_chains[i] + '__embed"/>\n')
		f.write('            <Data desc="' + mp_chains[i] + '_span"  resource_tag="' + mp_chains[i] + '__span"/>\n')
		f.write('            <Data desc="' + mp_chains[i] + '" resource_tag="' + mp_chains[i] + '_"/>\n')

	# Finish up the jobs and the file
	f.write('        </Job>\n')
	f.write('    </Jobs>\n')
	f.write('</JD2ResourceManagerJobInputter>\n')

	# Close the resulting file
	f.close()

def main( argv ):

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


