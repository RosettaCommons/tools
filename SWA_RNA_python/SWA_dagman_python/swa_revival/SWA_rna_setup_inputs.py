#!/usr/bin/env python

import urllib
import subprocess
import argparse
from os.path import exists

###################################################################################################

def fetch_pdb(id):
	url = 'http://www.rcsb.org/pdb/files/%s.pdb' % id
	return urllib.urlopen(url).read()

def download_pdb(id):
	pdb=fetch_pdb(id.upper())
	if pdb:
		pdb_out=id+'.pdb'
		pdb_fout=open(pdb_out,'w')
		pdb_fout.write(pdb)
		pdb_fout.close()
		return pdb_out
	else:
		return None

def get_scaffold_tag(pdb_file, sample_res):
	nts=''
	for sr in sample_res:
		for line in open(pdb_file, 'r').readlines():
			if line.split()[5] == sr:
				nts+=line.split()[3]
				break	
	scaffold_tag='no'+nts+'_'
	return scaffold_tag

def extract_chain(pdb_file, chain_id):
	pdblines=[line for line in open(pdb_file, 'r').readlines() if ('ATOM' in line or 'TER' in line) and chain_id.upper() in line] 
	pdb_file_fixed=pdb_file.replace('.pdb',str('_'+chain_id+'.pdb'))
	pdb_fout=open(pdb_file_fixed,'w')
	for line in pdblines:
		pdb_fout.write(line)
	pdb_fout.close()
	return pdb_file_fixed

def make_rna_rosetta_ready(pdb_file):
	rna_pdb_in=pdb_file.replace('.pdb','_RNA.pdb')
	subprocess.call('make_rna_rosetta_ready.py '+pdb_file, shell=True)
	return rna_pdb_in

def setup_fasta(pdb_file):
	fasta_out=pdb_file.replace('.pdb','.fasta')
	subprocess.call('pdb2fasta.py '+pdb_file+' > '+fasta_out, shell=True)
	return

def setup_scaffold(pdb_file, sample_res):
	scaffold_tag=get_scaffold_tag(pdb_file, sample_res)
	subprocess.call('pdbslice.py '+pdb_file+' -excise '+' '.join(sample_res)+' '+scaffold_tag, shell=True)
	return


###################################################################################################3

#SWA_rna_setup_inputs.py <pdb_id> <sample_res_list> [-chain_id <chain_id>]
parser = argparse.ArgumentParser(description='')
parser.add_argument('pdb_id', help='Name of a pdb')
parser.add_argument('-sample_res', nargs='+', help='Sample residues', default=None)
parser.add_argument('-chain_id', help='Chain ID', default=None)
parser.add_argument('-new_dir', default=False)
args=parser.parse_args()

pdb_id=(args.pdb_id).lower()
chain_id=(args.chain_id).lower()
sample_res=args.sample_res
new_dir=args.new_dir

###################################################################################################

print 'pdb_id = '+pdb_id
print 'chain_id = '+chain_id
print 'sample_res = ',sample_res

if new_dir:		
	subprocess.call('mkdir inputs', shell=True)
pdb_file=download_pdb(pdb_id)
if chain_id:	pdb_file=extract_chain(pdb_file, chain_id)
rna_pdb_in=make_rna_rosetta_ready(pdb_file)
setup_fasta(rna_pdb_in)
if sample_res:	setup_scaffold(rna_pdb_in, sample_res)


if exists('./swa/clean_directory'):
	subprocess.call('cp *'+pdb_id+'* '+'./swa/clean_directory/', shell=True )
if exists('./swm/clean_directory'):
	subprocess.call('cp *'+pdb_id+'* '+'./swm/clean_directory/', shell=True )
subprocess.call('rm ./*'+pdb_id+'*', shell=True)