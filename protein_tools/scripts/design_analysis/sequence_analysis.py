#!/usr/bin/env python
import sys
import math
import warnings
import os
import copy
try:
	path = "/blue/meilerlab/apps/scripts/libraries/python_library/"
	sys.path.append(path)
	import Bio.SeqIO as seqIO
	from Bio.PDB.PDBParser import PDBParser as pdbparse
	from Bio.PDB import PDBExceptions
except ImportError:
	paths = "    ".join(sys.path)
	raise ImportError(
		"We have tried finding a Biopython path for you, but it's just not in the cards, here is your search path:"+paths)

def three_to_one(amino_acid):
	amino_acid_names = { "ALA" : "A","GLY":"G","LEU":"L","ILE":"I","MET":"M","PHE":"F","VAL":"V","SER":"S","PRO":"P","THR":"T","TYR":"Y","HIS":"H","GLU":"E",\
						"ASP": "D","ASN":"N","GLN":"Q","TRP":"W","CYS":"C","CYD":"C","ARG":"R","LYS":"K","TES":"X","S56":"X"}
        if amino_acid not in amino_acid_names:
            return 'X'
        else:
	    return amino_acid_names[amino_acid]

# prelim junk#
warnings.simplefilter('ignore', PDBExceptions.PDBConstructionWarning)
# end prelim junk#

class DesignParsers:

	def parse_fasta_single(self, fasta):
		"""return a biopython fasta object"""
		try:
			return seqIO.read(open(fasta), "fasta")
		except IOError:
			raise IOError("can't find "+str(fasta)+" it may not exist in this path")

	def parse_pdbs(self, pdbs):
		"""Return a a list of biopython pdb objects\n doesn't matter if you pass several or only one pdb"""
		p = pdbparse(PERMISSIVE=1)
		try:
			if type(pdbs) == str:
				return p.get_structure(pdbs, pdbs)
			else:
				return [p.get_structure(i, i) for i in pdbs]
		except:
			raise IOError("can't find "+str(pdbs)+" it may not exist in this path")

	def parse_corr(self, corr):
		corr_dict = {}
		try:
			with open(corr) as f:
				for line in f.readlines():
					entity = line.split()[0]
					res_number = line.split()[1]
					chain = line.split()[2]
					corr_dict[entity] = {"chain": chain,"res":res_number}
			return corr_dict
		except IOError:
			print "corr file didn't parse"
			return ""

	def parse_res(self, res):
		self.res_dict = {}
		try:
			with open(res) as f:
				breaker = False
				for z, line in enumerate(f.readlines()):
					if line.split()[0] == "start":
						breaker = True
						continue
					if breaker:
						entity = z
						res_number = line.split()[0]
						chain = line.split()[1]
						self.res_dict[entity] = {"chain": chain,"res":res_number}
			return self.res_dict
		except IOError:
			print "res file didn't parse"
			return ""


class Customdict(dict):
	pass


class DesignScoreTable():

	"""A class encasulating the entire design scoring of every pdb. The __init__ does the initial analysis depending on what options you provide"""

	def __init__(self,native_pdb="",native_fasta="",t_type="",targets=[],corr_file="",rez_file="",out_files=[]):
			# arguments we will encapsulate
		self.t_type = t_type
		self.corr_file = corr_file
		self.res_file = rez_file
		self.native_pdb = native_pdb
		self.native_fasta = native_fasta
		self.targets = targets
		self.designed_keys = []
		self.native_keys = {}

		# analysis dictionaries we will only fill them if the options specify
		self.analysis_dict = {}
			# we are housing a whole lot of analysis dicts here, their primary keys
			# are always a tuple ("chain","res_number") analysis_dict will compare
			# every pdb of fasta file and get the occurences of everything
		self.designed_dict = {}
			#designed_dict  will be a subset of analysis_dict but only for the residues we defined in a res or corr file
		self.native_dict = {}
			#native_dict will be a subdict of analysis dict but only showing how many times the native amino acid was designed
		self.native_dict_designed = {
			}  # no doubt the most important dict, a subset of designed_dict_bit but only with native and amino acid bit scores shown

		#Run methods
		if 'a' in out_files:
			self._analyze_all(self.t_type,self.targets)

		if 'd' in out_files:
			if self.corr_file:
				self._get_designed_keys(corr_file=self.corr_file)
			if self.res_file:
				self._get_designed_keys(corr_file=self.res_file)
			self._analyze_all(self.t_type, self.targets)

		if 'n' in out_files or 'nd' in out_files:
			if self.corr_file:
				self._get_designed_keys(corr_file=self.corr_file)
			if self.res_file:
				self._get_designed_keys(corr_file=self.res_file)
			self._get_native_pdb_keys()
			self._analyze_all(self.t_type, self.targets)

		if 'nd' in out_files:
			self._get_native_designed_residues()

	# BEGIN MEMBER HELPER FUNCTIONS####
	def _analyze_all(self, target_type, targets):
		"""The primary method that gets the occurence matrix and fills analysis_dict"""
		de = DesignScoreEntry()  # calls on the design score entry class which will analyze one pdb at a time
		for i in targets:
			current_dictionary = de.define_target(
				target_type, i)  # define all the pdbs into a member variable called target_pdbs, contains file names too
		occurences = de.get_occurences(
			current_dictionary)  # get the occurences back, basically mash all the tables together and ask how often we see a particular amino acid

		#used as a holder so we don't initialize the dictionary twice
		for key in occurences:
			chain = key[0]  # break down the occurence dictionary so we can actually explain it
			rez_number = key[1]
			aa = key[2]  # what is the aa and how many times does it occur
			occur = occurences[key]
			key = (chain, rez_number)  # the key will will use


			#if aa is found in occurences, assign it. Or else don't assign it. We will catch the exception later on
			if self.native_keys:
				if self.native_keys[key] == aa:  # see if the native_key amino acid equals the amino acid and that its not
					self.native_dict[key] = (self.native_keys[key], occur)


			if self.designed_keys:
				if key in self.designed_keys or key[0] == "null":
					if key not in self.designed_dict:  # see if that key is in their, if not we fill it with all the amino acids with 0 occurences
						self.designed_dict[key] = {'A': 0,'V':0,'G':0,'L':0,'I':0,'W':0,'P':0,'C':0,'N':0,'Q':0,'E':0,'D':0,'Y':0,'R':0,'H':0,'M':0,'S':0,'T':0,'K':0,'F':0} #fill dictionary with zero occurences of all amino acids
					self.designed_dict[key][aa] = int(
						occur)  # fill each chain, and rez number with the amount of times that aa was present

			else:
				if key not in self.analysis_dict:  # see if that key is in their, if not we fill it with all the amino acids with 0 occurences
					self.analysis_dict[key] = {'A': 0,'V':0,'G':0,'L':0,'I':0,'W':0,'P':0,'C':0,'N':0,'Q':0,'E':0,'D':0,'Y':0,'R':0,'H':0,'M':0,'S':0,'T':0,'K':0,'F':0} #fill dictionary with zero occurences of all amino acids
				self.analysis_dict[key][aa] = int(
					occur)  # fill each chain, and rez number with the amount of times that aa was present

	def _get_designed_keys(self,corr_file="",res_file=""):
		if corr_file:
			for i in corr_file:
				chain = corr_file[i]['chain']
				rez = int(corr_file[i]['res'])
				key = (chain, rez)
				self.designed_keys.append(key)

	def get_bit_scores(self,normalize=False,dict=""):
		"""This method is extremely useful as it will turn the numbers in analysis dict into bitscores, option to normalize or not with normalized=True)"""
		_return_dict = {}
		for line in dict:  # each line is a residue
			sum_l = 0.0  # or total number of amino acids at any given positon, this is used for the propensity calculation
			p_entropy = 0.0  # entropy calculates how much a position verfied
			sum_l = sum([float(x) for x in dict[line].values()])  # add up all the amino acid occurences, should equal the number of models or fastas you passed
			_return_dict[line] = dict[line]  # copy the key from the parameter dict
			for aa in dict[line]:  # iterate through the amino acids
				bit_score = 0.0
				try:
					propensity = float(dict[line][aa])/float(sum_l)  # propensity is how many times you see that particular amino acid / how many amino acids total at that position
					bit_score = float(propensity)*math.log((20.00*float(
						propensity)), 2.0)  # propensity * log2(propensity*how many amino acids where considered)
				except ValueError:
					bit_score = 0.0  # if we never saw an amino acid its bit score is zero
				_return_dict[line][aa] = float(("%.2f") % bit_score)  # bit score dictionary
			if normalize:
				p_entropy = sum(_return_dict[line].values())
				for aa in _return_dict[line]:
					_return_dict[line][aa] = ("%.2f") % ((_return_dict[line][aa]/p_entropy))  # normalize for each position by dividing it by the entropy (that is the maximum acheivable score)
			else:
				for aa in _return_dict[line]:
					# print aa, _return_dict[line][aa]
					_return_dict[line][aa] = ("%.2f") % (_return_dict[line][aa])  # else just give back the regular bit scores
		return _return_dict

	def _get_native_fasta_keys(self):

		chain = "null"
		for i, j in enumerate(self.native_fasta.seq, start=1):
			self.native_keys[(chain, i)] = j

	def _get_native_pdb_keys(self):

		for i in self.native_pdb.get_residues():
			self.native_keys[(i.get_parent().get_id(), int(i.get_id()[1]))] = three_to_one(i.resname)

	def _get_native_designed_residues(self):
		"""pulls native residues defined in a score file out of another dictionary"""
		for i in self.designed_keys:
			#this is such garbage, but if there is an occurence of our native amino acid, we assign the occurence number, if the key isn't found, we never found that amino acid. This comes from the analyze all method
			try:
				self.native_dict_designed[i] = (self.native_keys[i], self.native_dict[i][1])
			except KeyError:
				self.native_dict_designed[i] = (self.native_keys[i],0)
	def get_native_designed_residues(self,dict=""):
		"""pulls native residues defined in a score file out of another dictionary"""
		_return_dict = {}
		for i in self.designed_keys:
			_return_dict[i] = (self.native_keys[i], dict[i][self.native_keys[i]])
		return _return_dict

	def get_native_residues(self,dict=""):
		"""pulls native residues out of another dictionary"""
		_return_dict = {}
		for i in dict:
				_return_dict[i] = (self.native_keys[i], dict[i][self.native_keys[i]])
		return _return_dict

	def _verbosely_print_dict(self, dict):
		"""my super sketchy method of printing out this analysis dictionary in order"""
		header = list(max(dict.values()))
		print '    ',
		for d in sorted(header):
			print d.center(3),
		print
		for key in sorted(dict):
			print key[0]+","+str(key[1])+' ',
			for y in sorted(dict[key]):
				print str(dict[key][y]).center(3),
			print

	def _verbosely_print_native_dict(self, dict):
		for line in sorted(dict):
			print "chain {:<1},res_number {:>2},{:^1} {:<1} {:<1} {:<1}".format(line[0], line[1], "native", dict[line][1], "===>", dict[line][0])

	def _write_out_dict(self,filename="",dict="",prefix="",length_of_pdbs=0):
		space = "tab"
		with open(prefix+filename+"."+space, 'w') as f:
			header = list(max(dict.values()))
			f.write("Chain,Position\t")
			for entry in sorted(header):
				string = "{:<3}".format(entry)
				f.write(string)
			f.write("\n")
			for key in sorted(dict):
				chain = key[0]
				rez = key[1]
				string = chain + "," + str(rez)
				f.write("{:^15}".format(string))
				for d in sorted(dict[key]):
					if length_of_pdbs:
						occur = float(dict[key][d])/length_of_pdbs
					else:
						occur = float(dict[key][d])
					f.write("{:^5.3f} ".format(occur).center(5))
				f.write("\n")

	def _write_out_table(self,filename="",dict="",prefix="",space="",length_of_pdbs=0):
		space = "tab"
		with open(prefix+filename+"."+space, 'w') as f:
			f.write("{0}\t{1}\t{2}\n".format("chain,position", "occurences", "native amino acid"))
			for line in sorted(dict):
				occur = float(dict[line][1])/length_of_pdbs
				chain = line[0]
				residue = line[1]
				identity = dict[line][0]
				#print occur,chain,residue,identity
				f.write("{0},{1}\t".format(line[0],line[1])+"{0:.2f}\t{1}\n".format(occur,dict[line][0]))


class DesignScoreEntry():

	"""A class Encapsulating the scores of a particular pdb"""
	def __init__(self):
		self.target_dict = {}  # contains all of the amino acids in question

	def define_target(self, target_type, name):
		if target_type == 'fasta':
			self.fasta_name = name.id.split('.fasta')[0]
			for i, g in enumerate(name.seq, start=1):
				self.target_dict[self.fasta_name, "null", i] = g
		else:
			self.pdb_name = os.path.basename(name.get_id())
			for res in name.get_residues():
				self.chain_name = res.get_full_id()[
												  2]  # IMPORTANT residue.get_full_id returns tuple ("pdb name,modelnum,chain,('',residue_num,''))
				self.residue_num = res.get_full_id()[3][
												   1]  # IMPORTANT residue.get_full_id returns tuple ("pdb name,modelnum,chain,('',residue_num,''))
				self.target_index = (self.pdb_name, self.chain_name, self.residue_num)  # the index is our absolute runique key
				self.target_dict[self.target_index] = three_to_one(res.resname)
		return self.target_dict

	def define_fasta(self, fasta_name):
		self.fasta_name = fasta_name.id.split('.fasta')[0]
		for i, g in enumerate(fasta_name.seq, start=1):
			self.target_dict[self.fasta_name, "null", i] = g

	def define_native(self, pdb_name):
		self.native_dict = self.define_target(pdb_name)

	def get_occurences(self, target_dict):
		"""where the magic happens"""
		occurence_dict = {}
		for i in target_dict:
			chain = i[1]
			seq = i[2]
			residue = target_dict[i]
			try:
				occurence_dict[chain, seq, residue] += 1
			except KeyError:
				occurence_dict[chain, seq, residue] = 1
		return occurence_dict

	def get_fasta_occurences(self):
		"""where the magic happens"""
		occurence_dict = {}
		for i in self.target_dict:
			file = i[0]
			resnum = i[1]
			chain = i[2]
			residue = self.target_dict[i]
			try:
				occurence_dict[chain, resnum, residue] += 1
			except KeyError:
				occurence_dict[chain, resnum, residue] = 1
		return occurence_dict
