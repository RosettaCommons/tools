#!/usr/bin/python

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
		subprocess.call('rm '+file, shell=True)
	return


###################################################################################################


class PDBInputFactory(object):


	def init(self, pdb_id, chain_id, sample_res_str):
		self.pdb_id_=str(pdb_id).lower()
		self.chain_id_=str(chain_id).lower()
		self.sample_res_str_=sample_res_str
				
		self.native_pdb=None
		self.fasta_file=None
		self.scaffold_file=None


	def setup_pdb_inputs(self, pdb_id, chain_id, sample_res_str):
		self.init(pdb_id, chain_id, sample_res_str)
		self.download_pdb()
		self.extract_chain()
		self.setup_fasta()
		self.setup_scaffold()
		self.make_rna_rosetta_ready()
		#slice_sample_res_and_surrounding(self):
		return


	def fetch_pdb(self):
		url = 'http://www.rcsb.org/pdb/files/%s.pdb' % self.pdb_id_.upper()
		print 'Getting ... '+url 
		return urllib.urlopen(url).read()


	def download_pdb(self):
		self.native_pdb=self.pdb_id_+'.pdb'
		open(self.native_pdb,'w').write(self.fetch_pdb())
		print 'Writing ... '+self.native_pdb
		return


	def extract_chain(self):
		if self.chain_id_:
			native_pdb_fixed=self.native_pdb.replace('.pdb',str('_'+self.chain_id_+'.pdb'))
			pdblines=''
			for line in open(self.native_pdb, 'r').readlines():
				if( (line.split()[0]=='ATOM' and line.split()[4]==self.chain_id_.upper()) or 
					(line.split()[0]=='TER'  and line.split()[3]==self.chain_id_.upper()) ):
					pdblines+=line 
			open(native_pdb_fixed,'w').write(pdblines)
			self.native_pdb=native_pdb_fixed
			print 'Writing ... '+self.native_pdb
		return


	def setup_fasta(self):
		self.fasta_file=self.native_pdb.replace('.pdb','.fasta')
		print "Running ... 'pdb2fasta.py "+self.native_pdb+" > "+self.fasta_file+"'"
		subprocess.call(   'pdb2fasta.py '+self.native_pdb+' > '+self.fasta_file, shell=True  )
		return


	def setup_scaffold(self):
		nts=''
		for sr in self.sample_res_str_.split(' '):
			for line in open(self.native_pdb, 'r').readlines():
				if line.split()[5] == sr:
					nts+=line.split()[3]
					break	
		scaffold_tag='no'+nts+'_'
		self.scaffold_file=scaffold_tag+self.native_pdb
		print "Running ... 'pdbslice.py "+self.native_pdb+" -excise "+self.sample_res_str_+" "+scaffold_tag+"'"
		subprocess.call(   'pdbslice.py '+self.native_pdb+' -excise '+self.sample_res_str_+' '+scaffold_tag, shell=True  )
		return


	def make_rna_rosetta_ready(self):
		print "Running ... 'make_rna_rosetta_ready.py "+self.native_pdb+"'"
		subprocess.call(   'make_rna_rosetta_ready.py '+self.native_pdb, shell=True  )
		self.native_pdb=self.native_pdb.replace('.pdb','_RNA.pdb')
		return

	def slice_sample_res_and_surrounding(self):
		command='~/src/rosetta/tools/SWA_RNA_python/SWA_dagman_python/misc/slice_sample_res_and_surrounding.py'
		command+='-loop_segment_list %s'%str(self.sample_res_str_.split(' ')[0]+'-'+self.sample_res_str_.split(' ')[-1])
		command+='-pose_name %s'%self.native_pdb
		command+='-MODE manual_sub'
		subprocess.call( command, shell=True  )
		return

###################################################################################################


class JobFileFactory(object):

	def __init__(self, use_oldscore=True):
		### SWM README FILE
		self.swm_out_file_silent='swm_rebuild_oldscore.out'
		self.nstruct=10
		self.ncycles=200
		### SWM SUBMIT FILE
		self.swm_out_dir='oldscore_out'
		self.swm_num_slave_nodes=50
		self.swm_master_wall_time=72
		self.use_oldscore=use_oldscore
		### SWA SUBMIT File
		self.swa_num_slave_nodes=300 
		self.swa_master_wall_time=72
		self.swa_master_memory_reserve=2048 
		self.swa_dagman_file='rna_build.dag'

		### scorefuncion specifics for SWM
		self.use_oldscore=use_oldscore
		if not self.use_oldscore:
			self.swm_out_dir='newscore_out'
			self.swm_out_file_silent='swn_rebuild_newscore.out'


	def init(self, directory, native_pdb, fasta_file, scaffold_file, sample_res_str ):
		self.directory=directory
		self.native_pdb=native_pdb
		self.fasta_file=fasta_file
		self.scaffold_file=scaffold_file
		self.sample_res_str=sample_res_str																            #40 41 42 43 44 45
		self.all_res_str=str( self.get_first_res()+'-'+self.get_last_res() )   												#1-77 
		self.input_res_str=str( self.get_first_res()+'-'+str(int(sample_res_str.split(' ')[0])-1)+' '
								+str(int(sample_res_str.split(' ')[-1])+1)+'-'+self.get_last_res()	) 	#1-39 45-77
		

	def get_first_res(self):
		for line in open(self.directory+'/'+self.scaffold_file, 'r').readlines():
			if line.split()[0]=='ATOM':
				first_res_idx=line.split()[5]
				break
		return first_res_idx

	def get_last_res(self):
		for line in reversed(open(self.directory+'/'+self.scaffold_file, 'r').readlines()):
			if line.split()[0]=='ATOM':
				last_res_idx=line.split()[5]
				break		
		return last_res_idx

	def setup_clean_directory(self, directory, native_pdb, fasta_file, scaffold_file, sample_res_str ):
		self.init(directory, native_pdb, fasta_file, scaffold_file, sample_res_str )
		self.write_README(directory)
		self.write_SUBMIT(directory)


	def write_README(self, directory):
		command=None
		if 'swm' in directory:
			command='~/src/rosetta/main/source/bin/stepwise \\\n' 
			command+='-out:file:silent %s \\\n'%self.swm_out_file_silent
			if self.use_oldscore:
				command+='-score:rna_torsion_potential ps_03242010/ \\\n'
				command+='-score:weights stepwise/rna/rna_loop_hires_04092010.wts \\\n'
				command+='-analytic_etable_evaluation false \\\n'
			command+='-constant_seed \\\n'
			command+='-nstruct %d \\\n'%self.nstruct  
			command+='-cycles %d \\\n'%self.ncycles
			command+='-native %s \\\n'%self.native_pdb
			command+='-fasta %s \\\n'%self.fasta_file 
			command+='-s %s \\\n'%self.scaffold_file
			command+='-sample_res %s \\\n'%self.sample_res_str
			command+='-rmsd_res %s \\\n'%self.sample_res_str
			command+='-global_sample_res_list %s \\\n'%self.sample_res_str  
			command+='-cutpoint_closed %s \\\n'%self.sample_res_str
			command+='-input_res %s \\\n'%self.input_res_str		#1-39 45-77
			command+='-stepwise:fixed_res %s \\\n'%self.input_res_str 
			command+='-alignment_res %s \\\n'%self.input_res_str
			command+='-jump_point_pairs %s \\\n'%self.all_res_str 
		if 'swa' in directory:
			command='~/src/rosetta/tools/SWA_RNA_python/SWA_dagman_python/SWA_DAG/setup_SWA_RNA_dag_job_files.py \\\n'
			command+='-single_stranded_loop_mode True \\\n'  
			command+='-native_pdb %s \\\n'%self.native_pdb					#3d2v_b_RNA.pdb'
			command+='-fasta %s \\\n'%self.fasta_file						#3d2v_b_RNA.fasta' 
			command+='-s %s \\\n'%self.scaffold_file							#noUUGAA_3d2v_b_RNA.pdb'  
			command+='-sample_res %s \\\n'%self.sample_res_str					#40 41 42 43 44' 
		if not command:	return None 
		readme_file=directory+'/README'
		print 'Writing ... '+readme_file
		open(readme_file,'w').write(command) 
		return readme_file


	def write_SUBMIT(self, directory):
		command=None
		if 'swm' in directory:
			command='rosetta_submit.py \\\n'
			command+='./README \\\n'
			command+='%s \\\n'%self.swm_out_dir
			command+='%d \\\n'%self.swm_num_slave_nodes
			command+='%d \\\n'%self.swm_master_wall_time
		if 'swa' in directory:
			command='~/src/rosetta/tools/SWA_RNA_python/SWA_dagman_python/dagman/submit_DAG_job.py \\\n' 
			command+='-master_wall_time %d \\\n'%self.swa_master_wall_time 
			command+='-master_memory_reserve %d \\\n'%self.swa_master_memory_reserve 
			command+='-num_slave_nodes %d \\\n'%self.swa_num_slave_nodes 
			command+='-dagman_file %s \\\n'%self.swa_dagman_file
		if not command:	return None 
		submit_file=directory+'/SUBMIT'
		print 'Writing ... '+submit_file
		open(submit_file, 'w').write(command)
		return submit_file


###################################################################################################