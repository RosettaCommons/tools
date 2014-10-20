#!/usr/bin/python

import urllib
import subprocess
import argparse
from os.path import exists

###################################################################################################

class PDBInputFactory(object):

	def __init__(self):
		return

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
		subprocess.call('make_rna_rosetta_ready.py '+self.pdb_file, shell=True)
		return rna_pdb_file

	def setup_fasta(self):
		if self.is_rna:	pdb_file=self.rna_pdb_file
		else:	pdb_file=self.pdb_file
		fasta_out=pdb_file.replace('.pdb','.fasta')
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
		subprocess.call('pdbslice.py '+pdb_file+' -excise '+' '.join(self.sample_res_list_)+' '+scaffold_tag, shell=True)
		return scaffold_file


