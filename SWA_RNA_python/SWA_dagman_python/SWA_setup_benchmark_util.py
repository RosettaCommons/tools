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
		if exists(file):
			subprocess.call('rm '+file, shell=True)
	return


###################################################################################################


class PDBInputFactory(object):


	def init(self, pdb_id, chain_id, sample_res_list, is_rna):
		self.pdb_id_=str(pdb_id).lower()
		if chain_id:
			self.chain_id_=str(chain_id).lower()
		else:
			self.chain_id_=None
		if sample_res_list:
			self.sample_res_list_=sample_res_list
		else:
			self.sample_res_list_=None
		self.is_rna=is_rna
		
		self.pdb_file=None
		self.pdb_file_original=None
		self.rna_pdb_file=None
		self.fasta_file=None
		self.scaffold_file=None


	def setup_pdb_inputs(self, pdb_id, chain_id=None, sample_res_list=None, is_rna=True):
		self.init(pdb_id, chain_id, sample_res_list, is_rna)
		self.pdb_file=self.download_pdb()
		if self.pdb_file:
			self.pdb_file=self.extract_chain()
			if self.is_rna:
				self.rna_pdb_file=self.make_rna_rosetta_ready()
			self.fasta_file=self.setup_fasta()
			self.scaffold_file=self.setup_scaffold()
		return self.pdb_file_original, self.pdb_file, self.rna_pdb_file, self.fasta_file, self.scaffold_file


	def fetch_pdb(self):
		print 'Fetching ... '+self.pdb_id_.upper()+'.pdb'
		url = 'http://www.rcsb.org/pdb/files/%s.pdb' % self.pdb_id_.upper()
		return urllib.urlopen(url).read()


	def download_pdb(self):
		pdb=self.fetch_pdb()
		if pdb:
			pdb_file=self.pdb_id_+'.pdb'
			pdb_fout=open(pdb_file,'w')
			pdb_fout.write(pdb)
			pdb_fout.close()
			self.pdb_file_original=pdb_file
			return pdb_file
		else:
			return None


	def extract_chain(self):
		if self.chain_id_:
			print 'Extracting ... chain '+self.chain_id_.upper()
			pdb_file_fixed=self.pdb_file.replace('.pdb',str('_'+self.chain_id_+'.pdb'))
			pdblines=[line for line in open(self.pdb_file, 'r').readlines() if ('ATOM' in line or 'TER' in line) and self.chain_id_.upper() in line] 
			pdb_fout=open(pdb_file_fixed,'w')
			for line in pdblines:	pdb_fout.write(line)
			pdb_fout.close()
			return pdb_file_fixed
		else:
			return self.pdb_file


	def make_rna_rosetta_ready(self):
		rna_pdb_file=self.pdb_file.replace('.pdb','_RNA.pdb')
		print 'Running ... make_rna_rosetta_ready.py '+self.pdb_file
		subprocess.call('make_rna_rosetta_ready.py '+self.pdb_file, shell=True)
		return rna_pdb_file


	def setup_fasta(self):
		if self.is_rna:	pdb_file=self.rna_pdb_file
		else:	pdb_file=self.pdb_file
		fasta_out=pdb_file.replace('.pdb','.fasta')
		print 'Running ... pdb2fasta.py '+pdb_file+' > '+fasta_out
		subprocess.call('pdb2fasta.py '+pdb_file+' > '+fasta_out, shell=True)
		return fasta_out


	def setup_scaffold(self):
		if self.is_rna:	pdb_file=self.rna_pdb_file
		else:	pdb_file=self.pdb_file
		nts=''
		for sr in self.sample_res_list_:
			for line in open(pdb_file, 'r').readlines():
				if line.split()[5] == sr:
					nts+=line.split()[3]
					break	
		scaffold_tag='no'+nts+'_'
		scaffold_file=scaffold_tag+pdb_file
		print 'Running ... pdbslice.py '+pdb_file+' -excise '+' '.join(self.sample_res_list_)+' '+scaffold_tag
		subprocess.call('pdbslice.py '+pdb_file+' -excise '+' '.join(self.sample_res_list_)+' '+scaffold_tag, shell=True)
		return scaffold_file


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
		self.native_pdb=native_pdb
		self.fasta_file=fasta_file
		self.scaffold_file=scaffold_file
		self.sample_res_str=sample_res_str 																	#40 41 42 43 44 45
		self.all_res_str=str(open(directory+'/'+fasta_file, 'r').readline().split(':')[-1].replace('\n',''))			#1-77 from fasta
		self.input_res_str=str(  self.all_res_str.split('-')[0]+'-'+str(int(sample_res_str.split(' ')[0])-1)+' '
								+str(int(sample_res_str.split(' ')[-1])+1)+'-'+self.all_res_str.split('-')[-1]		) 	#1-39 45-77
		

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