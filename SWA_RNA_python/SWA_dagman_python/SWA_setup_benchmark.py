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
def setup_clean_directory(directory):
	write_README(directory)
	write_flags(directory)
	write_SUBMIT(directory)

def write_README(directory):
	command=''
	if 'swm' in directory:
		command+='~/src/rosetta/main/source/bin/stepwise' 
		command+=' @flags' 
		command+=' -out:file:silent swm_rebuild_oldscore.out'
		command+=' -score:rna_torsion_potential ps_03242010/'
		command+=' -score:weights stepwise/rna/rna_loop_hires_04092010.wts'
		command+=' -analytic_etable_evaluation false'
	elif 'swa' in directory:
		command+='~/src/rosetta/tools/SWA_RNA_python/SWA_dagman_python/SWA_DAG/setup_SWA_RNA_dag_job_files.py'
		command+=' @flags'
	else:
		return None	
	readme_file=directory+'/README'
	readme=open(readme_file,'w')
	readme.write(command)
	readme.close()
	return readme_file 

def write_flags(directory):
	flags=[]
	if 'swm' in directory:
		flags+='-s %s\n'%scaffold_file
		flags+='-native %s\n'%native_pdb
		flags+='-fasta %s\n'%fasta_file 
		flags+='-nstruct %d\n'%nstruct  
		flags+='-constant_seed\n'
		flags+='-cycles %d\n'%ncycles
		flags+='-stepwise:fixed_res  %d-%d %d-%d\n'%fixed_res_ranges #1-39 45-77 
		flags+='-rmsd_res  40 41 42 43 44\n' 
		flags+='-jump_point_pairs  1-77\n'  
		flags+='-alignment_res  1-39 45-77\n' 
		flags+='-global_sample_res_list 40 41 42 43 44\n'    
		flags+='-sample_res  40 41 42 43 44\n'
		flags+='-cutpoint_closed 40 41 42 43 44\n' 
		flags+='-input_res 1-39 45-77\n'
	elif 'swa' in directory:
	
	else:
		return None	
	flags_file_name=directory+'/flags'
	flags_file=open(flags_file_name,'w')
	flags_file.write(flags)
	flags_file.close()
	return flags_file_name 

def write_SUBMIT(directory):
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


	# Set up README and SUBMIT files
	setup_clean_directory(swm_clean_dir)
	setup_clean_directory(swm_clean_dir)

###################################################################################################