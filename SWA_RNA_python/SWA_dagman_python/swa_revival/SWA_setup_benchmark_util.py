#!/usr/bin/env python

import urllib
import subprocess
import argparse
from os.path import exists

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
		subprocess.call('rm -r '+file, shell=True)
	return


###################################################################################################


class PDBInputFactory(object):


	def init(self, pdb_id, chain_id, sample_res_str):
		self.pdb_id_=str(pdb_id).lower()
		self.chain_id_=str(chain_id).lower()
		self.sample_res_str_=sample_res_str
				
		self.native_pdb=None
		self.fasta_file=None
		self.template_pdb=None

		self.surrounding_radius=5.0
		self.no_renumber=False
		self.verbose=True


	def setup_pdb_inputs(self, pdb_id, chain_id, sample_res_str):
		self.init(pdb_id, chain_id, sample_res_str)
		self.download_pdb()
		self.clean_rna()
		#self.make_rna_rosetta_ready()
		self.extract_chain()
		self.setup_fasta()
		self.slice_sample_res_and_surrounding()
		self.setup_template_pdb()
		

	def subprocess_call_command(self, command):
		if self.verbose:	print "Running ... '"+command+"'"
		subprocess.call( command, shell=True )
		

	def fetch_pdb(self):
		url = 'http://www.rcsb.org/pdb/files/%s.pdb' % self.pdb_id_.upper()
		if self.verbose:	print 'Getting ... '+url 
		return urllib.urlopen(url).read()


	def download_pdb(self):
		self.native_pdb=self.pdb_id_+'.pdb'
		open(self.native_pdb,'w').write(self.fetch_pdb())
		if self.verbose:	print 'Writing ... '+self.native_pdb

	
	def clean_rna(self):
		pdb_in=open(self.native_pdb,'r').readlines()
		open(self.native_pdb,'w').write(''.join([x for x in pdb_in if x.split()[0]=='ATOM' and x.split()[3] in 'AUGC'])) 
		if self.verbose:	print 'Cleaning ... '+self.native_pdb


	def extract_chain(self):
		pdb_in=open(self.native_pdb,'r').readlines()
		self.native_pdb=self.native_pdb.replace(self.pdb_id_, str(self.pdb_id_+'_'+self.chain_id_))
		open(self.native_pdb,'w').write(''.join([x for x in pdb_in if x.split()[0]=='ATOM' and x.split()[4][0]==self.chain_id_.upper()])) 
		if self.verbose:	print 'Writing ... '+self.native_pdb
		

	def setup_fasta(self):
		self.fasta_file=self.native_pdb.replace('.pdb','.fasta')
		command="pdb2fasta.py "+self.native_pdb+" > "+self.fasta_file
		self.subprocess_call_command(command)
		

	def get_template_pdb_tag(self):
		nts=''
		for sr in self.sample_res_str_.split(' '):
			for line in open(self.native_pdb, 'r').readlines():
				if( ( line.split()[5]==sr ) or 	   ### for <  3 digit res nums i.e., X 999
					( line.split()[4][1:]==sr ) ): ### for >= 4 digit res nums i.e., X1999
					nts+=line.split()[3]
					break	
		return 'no'+nts+'_'


	def setup_template_pdb(self):
		template_pdb_tag=self.get_template_pdb_tag()
		self.template_pdb=template_pdb_tag+self.native_pdb
		command='pdbslice.py '+self.native_pdb+' -excise '+self.sample_res_str_+' '+template_pdb_tag
		self.subprocess_call_command(command)
		

	def make_rna_rosetta_ready(self):
		command='misc/SWA_make_rna_rosetta_ready.py '+self.native_pdb
		self.native_pdb=self.native_pdb.replace('.pdb','_RNA.pdb')
		command+=' -output_pdb '+self.native_pdb
		if self.no_renumber: command+=' -no_renumber'
		self.subprocess_call_command(command)
		

	def slice_sample_res_and_surrounding(self):
		command='~/src/rosetta//main/source/bin/swa_rna_util.linuxgccrelease' 
		command+=' -algorithm slice_ellipsoid_envelope' 
		command+=' -database ~/src/rosetta//main/database/'
		command+=' -s %s'%self.native_pdb
		command+=' -sample_res %s'%self.sample_res_str_  
		command+=' -surrounding_radius %d'%self.surrounding_radius
		#tag='no_loop_ellipsoid_expand_radius_%d_'%self.surrounding_radius*10
		#self.native_pdb=tag+self.native_pdb
		#self.subprocess_call_command(command)
		open('SLICE_SAMPLE_RES_AND_SURROUNDING.src','w').write(command)


		

###################################################################################################


class JobFileFactory(object):

	def __init__(self, use_oldscore=True):
		### SWM README FILE
		self.swm_out_file_silent='swm_rebuild_oldscore.out'
		self.nstruct=5
		self.ncycles=200
		### SWM SUBMIT FILE
		self.swm_out_dir='oldscore_out'
		self.swm_num_slave_nodes=100
		self.swm_master_wall_time=8
		self.use_oldscore=use_oldscore
		### SWA SUBMIT File
		self.swa_num_slave_nodes=300 
		self.swa_master_wall_time=16
		self.swa_master_memory_reserve=2048 
		self.swa_dagman_file='rna_build.dag'

		### scorefuncion specifics for SWM
		self.use_oldscore=use_oldscore
		if not self.use_oldscore:
			self.swm_out_dir='newscore_out'
			self.swm_out_file_silent='swn_rebuild_newscore.out'


	def init(self, directory, native_pdb, fasta_file, template_pdb, sample_res_str ):
		self.directory=directory
		self.native_pdb=native_pdb
		self.fasta_file=fasta_file
		self.template_pdb=template_pdb
		self.sample_res_str=sample_res_str																            #40 41 42 43 44 45
		self.all_res_str=str( self.get_first_res()+'-'+self.get_last_res() )   												#1-77 
		self.input_res_str=str( self.get_first_res()+'-'+str(int(sample_res_str.split(' ')[0])-1)+' '
								+str(int(sample_res_str.split(' ')[-1])+1)+'-'+self.get_last_res()	) 	#1-39 45-77
		

	def get_first_res(self):
		for line in open(self.directory+'/'+self.template_pdb, 'r').readlines():
			if line.split()[0]=='ATOM':
				if len(line.split()[4]) > 1: 
					first_res_idx=line.split()[4][1:]
				else:
					first_res_idx=line.split()[5]
				break		
		return first_res_idx

	def get_last_res(self):
		for line in reversed(open(self.directory+'/'+self.template_pdb, 'r').readlines()):
			if line.split()[0]=='ATOM':
				if len(line.split()[4]) > 1: 
					last_res_idx=line.split()[4][1:]
				else:
					last_res_idx=line.split()[5]
				break		
		return str(last_res_idx)

	def setup_clean_directory(self, directory, native_pdb, fasta_file, template_pdb, sample_res_str ):
		self.init(directory, native_pdb, fasta_file, template_pdb, sample_res_str )
		self.write_README(directory)
		self.write_SUBMIT(directory)


	def write_README(self, directory):
		command=None
		if 'swm' in directory:
			command='~/src/rosetta/main/source/bin/stepwise' 
			command+=' -out:file:silent %s'%self.swm_out_file_silent
			if self.use_oldscore:
				command+=' -score:rna_torsion_potential ps_03242010/'
				command+=' -score:weights stepwise/rna/rna_loop_hires_04092010.wts'
				command+=' -analytic_etable_evaluation false'
			command+=' -constant_seed'
			command+=' -nstruct %d'%self.nstruct  
			command+=' -cycles %d'%self.ncycles
			command+=' -native %s'%self.native_pdb
			command+=' -fasta %s'%self.fasta_file 
			command+=' -s %s'%self.template_pdb
			command+=' -sample_res %s'%self.sample_res_str
			command+=' -rmsd_res %s'%self.sample_res_str
			command+=' -global_sample_res_list %s'%self.sample_res_str  
			command+=' -cutpoint_closed %s'%self.sample_res_str
			command+=' -input_res %s'%self.input_res_str		#1-39 45-77
			command+=' -stepwise:fixed_res %s'%self.input_res_str 
			command+=' -alignment_res %s'%self.input_res_str
			command+=' -jump_point_pairs %s'%self.all_res_str 
		if 'swa' in directory:
			command='~/src/rosetta/tools/SWA_RNA_python/SWA_dagman_python/SWA_DAG/setup_SWA_RNA_dag_job_files.py'
			command+=' -single_stranded_loop_mode True'  
			command+=' -native_pdb %s'%self.native_pdb					#3d2v_b_RNA.pdb'
			command+=' -fasta %s'%self.fasta_file						#3d2v_b_RNA.fasta' 
			command+=' -s %s'%self.template_pdb							#noUUGAA_3d2v_b_RNA.pdb'  
			command+=' -sample_res %s'%self.sample_res_str					#40 41 42 43 44' 
		if not command:	return None 
		readme_file=directory+'/README'
		print 'Writing ... '+readme_file
		open(readme_file,'w').write(command) 
		return readme_file


	def write_SUBMIT(self, directory):
		command=None
		if 'swm' in directory:
			command='rosetta_submit.py'
			command+=' ./README'
			command+=' %s'%self.swm_out_dir
			command+=' %d'%self.swm_num_slave_nodes
			command+=' %d'%self.swm_master_wall_time
		if 'swa' in directory:
			command='~/src/rosetta/tools/SWA_RNA_python/SWA_dagman_python/dagman/submit_DAG_job.py' 
			command+=' -master_wall_time %d'%self.swa_master_wall_time 
			command+=' -master_memory_reserve %d'%self.swa_master_memory_reserve 
			command+=' -num_slave_nodes %d'%self.swa_num_slave_nodes 
			command+=' -dagman_file %s'%self.swa_dagman_file
		if not command:	return None 
		submit_file=directory+'/SUBMIT'
		print 'Writing ... '+submit_file
		open(submit_file, 'w').write(command)
		return submit_file


###################################################################################################