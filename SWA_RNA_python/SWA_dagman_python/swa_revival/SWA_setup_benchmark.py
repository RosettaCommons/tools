#!/usr/bin/env python

from sys import exit
from os import getcwd
from os.path import exists
import subprocess
from SWA_setup_benchmark_util import *

###################################################################################################

targets_file='./targets'

try:
	exists(targets_file)
except:
	print 'ERROR: could not find '+targets_file
	exit(0)

TARGETS=[' '.join(line.split()) for line in open(targets_file, 'r').readlines()]
TARGET_DIRS=dict([(target,target.replace(' ','_')) for target in TARGETS])

###################################################################################################

pdbinputfactory=PDBInputFactory()
jobfilefactory=JobFileFactory()

###################################################################################################

for target in TARGETS:
	
	target_dir=TARGET_DIRS[target]
	safe_mkdir(target_dir)
	print '\nCreating benchmark target ... '+target_dir

	swm_dir=target_dir+'/swm'
	swa_dir=target_dir+'/swa'
	safe_mkdir(swm_dir)
	safe_mkdir(swa_dir)
	
	swm_clean_dir=swm_dir+'/clean_directory'
	swa_clean_dir=swa_dir+'/clean_directory'
	safe_mkdir(swm_clean_dir)
	safe_mkdir(swa_clean_dir)


	### setup PDB input files
	target_info=target.split()
	pdb_id=target_info[0]
	chain_id=target_info[1]
	sample_res_str=' '.join([ str(x) for x in xrange(int(target_info[2]),int(target_info[3])+1) ])
	if len(target_info) == 6:
		sample_res_str+=' '+' '.join([ str(x) for x in xrange(int(target_info[4]),int(target_info[5])+1) ])
	pdbinputfactory.setup_pdb_inputs( pdb_id, chain_id, sample_res_str )
	pdb_input_files=['./*.pdb','./*.fasta','./*.src']
	safe_cp(pdb_input_files, swm_clean_dir)
	safe_cp(pdb_input_files, swa_clean_dir)
	safe_rm(pdb_input_files)


	### setup README and SUBMIT files
	native_pdb=pdbinputfactory.native_pdb
	fasta_file=pdbinputfactory.fasta_file
	template_pdb=pdbinputfactory.template_pdb
	jobfilefactory.setup_clean_directory( swm_clean_dir, native_pdb, fasta_file, template_pdb, sample_res_str )
	jobfilefactory.setup_clean_directory( swa_clean_dir, native_pdb, fasta_file, template_pdb, sample_res_str )

###################################################################################################