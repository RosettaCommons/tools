#!/usr/bin/python

from sys import exit
from os import getcwd
from os.path import exists
import subprocess
from PDBInputFactory import *

###################################################################################################

def safe_mkdir(directory):
	if not exists(directory):
		subprocess.call('mkdir '+directory, shell=True)
	return

def safe_cp(files, directory):
	if exists(directory):
		for file in files:
			subprocess.call('cp '+file+' '+directory, shell=True)
	return

def safe_rm(files):
	for file in files:
		if exists(file):
			subprocess.call('rm '+file, shell=True)
	return

###################################################################################################

targets_file='./targets'

try:
	exists(targets_file)
except:
	print 'ERROR: could not find '+targets_file
	exit(0)

TARGETS=[' '.join(line.split()[:4]) for line in open(targets_file, 'r').readlines()]
TARGET_DIRS=dict([(target,target.replace(' ','_')) for target in TARGETS])

###################################################################################################

pdbinputfactory=PDBInputFactory()

###################################################################################################

for target in TARGETS:
	
	target_dir=TARGET_DIRS[target]
	safe_mkdir(target_dir)

	swm_dir=target_dir+'/swm'
	swa_dir=target_dir+'/swa'
	safe_mkdir(swm_dir)
	safe_mkdir(swa_dir)
	
	swm_clean_dir=swm_dir+'/clean_directory'
	swa_clean_dir=swa_dir+'/clean_directory'
	safe_mkdir(swm_clean_dir)
	safe_mkdir(swa_clean_dir)

	target_info=target.split()
	pdb_id=target_info[0]
	chain_id=target_info[1]
	sample_res_list=[ str(x) for x in xrange(int(target_info[2]),int(target_info[3])+1) ]

	pdb_input_files=pdbinputfactory.setup_pdb_inputs(	pdb_id, 
														chain_id=chain_id,
														sample_res_list=sample_res_list, 
														is_rna=True	 )

	safe_cp(pdb_input_files, swm_clean_dir)
	safe_cp(pdb_input_files, swa_clean_dir)

	safe_rm(pdb_input_files)

###################################################################################################