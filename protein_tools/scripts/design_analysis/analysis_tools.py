# tools for friendly access to the data initiated from sequence analysis

from operator import attrgetter
import copy
import sys
import pdb
try:
	# if on meilerlab enviroment
	path = "/blue/meilerlab/apps/scripts/libraries/python_library/"
	sys.path.append(path)
	import scoretable
	from sequence_analysis import three_to_one
except:
	paths = "\n".join(sys.path)
	raise ImportError(
		"rosettaScore_beta can not be found in this path\n" + paths)


class AnalysisTools():

	'''Analysis tools works hand-in-hand with sequence analysis, acessing some of the methoseq_analysis there.
	In fact, whenyou initialize it, it takes a sequence analysis class
	at=AnalysisTools(SequenceAnalysis(**kwarg) and gives you access to all of my analysis tools so far'''

	def __init__(self, seq_analysis, length_of_pdbs):
		'''Analysis tools takes in another class sequence analysis. It would be more elegant to superclass'''
		# takes in sequence analysis class object
		self.seq_analysis = seq_analysis
		self._arg_dict = {
			'a': self.seq_analysis.analysis_dict, 'd': self.seq_analysis.designed_dict,
			'n': self.seq_analysis.native_dict, 'nd': self.seq_analysis.native_dict_designed}
		self.len_of_args = length_of_pdbs

	def verbose(self, list_o):
		'''at.verbose returns nothing but prints to the screen whatever dictionary or list you called for anlysis on
				see score types on design_analysis.py --help for more'''
		for key in list_o:
			# arg dict has what the class calls each dictionary
			#if it a dictionary we print out a matrix of the dictionary
			if type(self._arg_dict[key].values()[0]) == dict:
				self.seq_analysis._verbosely_print_dict(self._arg_dict[key])
			# sometimes its not a dictinonary, just a list of native residues
			# with there occurences, there is another method for that
			else:
				# if it is a 1 line column just print out the column.
				self.seq_analysis._verbosely_print_native_dict(
				self._arg_dict[key])

	def output(self, prefix="", list_o="", bit="", normal=False):
		if bit is False:
			for key in list_o:
				if type(self._arg_dict[key].values()[0]) == dict:
					self.seq_analysis._write_out_dict(filename=key, prefix=prefix, dict=self._arg_dict[key], length_of_pdbs=self.len_of_args)
				else:
					self.seq_analysis._write_out_table(filename=key, prefix=prefix, dict=self._arg_dict[key],length_of_pdbs=self.len_of_args)
		else:
			for key in list_o:
				if type(self._arg_dict[key].values()[0]) == dict and key != 'n':
					bitscore_dict = self.seq_analysis.get_bit_scores(dict=self._arg_dict[key], normalize=normal)
					self.seq_analysis._write_out_dict(filename=key, prefix=prefix + "bit_", dict=bitscore_dict)
				else:
					if key == 'n':
						bitscore_dict = self.seq_analysis.get_bit_scores(dict=self.seq_analysis.analysis_dict, normalize=normal)
						bitscore_native = self.seq_analysis.get_native_residues(dict=bitscore_dict)
						self.seq_analysis._write_out_table(filename=key, prefix=prefix + "bit_", dict=bitscore_native,length_of_pdbs=1)
					if key == 'nd':
						bitscore_dict = self.seq_analysis.get_bit_scores(dict=self.seq_analysis.designed_dict, normalize=normal)
						bitscore_dict_native = self.seq_analysis.get_native_designed_residues(dict=bitscore_dict)
						self.seq_analysis._write_out_table(filename=key, prefix=prefix + "bit_", dict=bitscore_dict_native,length_of_pdbs=1)

	def get_stats(self, filename="output_stats.txt", bit=""):
		if bit:
			shannon_b = 0.0
			dicti = self.bitscore_dict
			for i in dicti:
				for b in dicti[i]:
					shannon_b += float(dicti[i][b])
			nbd = 0.0
			for i in self.bitscore_dict_native:
				nbd += float(self.bitscore_dict_native[i][1])
			print "Total Bit Score of Design ===> %.4f" % (nbd)
			print "Total Shannon Entropy of Design ===> %.4f" % (shannon_b)
			print "Normalized Bit Score for design ===> %.4f" % (nbd / shannon_b)
			with open(filename, 'w') as f:
				f.write("Total Bit Score of Design ===> {:.4}\nTotal Shannon Entropy of Design ===> {:.4}\nNormalized Bit Score for design ===> {:.4}".format(
					nbd, shannon_b, nbd / shannon_b))
		else:
			dicti = self.seq_analysis.native_dict_designed
			occur = 0.0
			total = 0.0
			counter = 0
			for line in sorted(dicti):
				occur = float(dicti[line][1])
				total += occur / self.len_of_args
				counter += 1
			recovery = (total / counter) * 100
			print "Sequence Recovery of Design ====> {:.4}%".format(recovery)
			with open(filename, 'w') as f:
				f.write(
					"Sequence Recovery of Design ====> {:.4}%".format(recovery))

			# f.write(line[0]+","+str(line[1])+"\t"+str(occur)+"\t"+str(dict[line][0])+"\n")
	def get_rosetta_energy_scores(self, d_of_interest, models):
		self.all_models = {}
		for i in models:
			rs = RosettaScores(self._arg_dict[d_of_interest], i)
			self.all_models[i] = rs.get_score()
		return self.all_models

	def print_rosetta_dict(self, args):
		print "{0},{1},{2},{3},{4}".format("model", "anlysis", "chain", "residue", "residue_number", "score")
		for p in dictionaries:
			d = get_rosetta_energy_scores(p, models, term)
			for i in sorted(d):
				for j in sorted(d[i]):
					printer = []
					for k in sorted(d[i][j]):
						printer.append(str(d[i][j][k]))
					print i.split(".")[0] + "," + ",".join(printer)

	def output_rosetta_dict(self, dictionaries, prefix, models, file="output"):
		name = prefix + "rosetta_scores.csv"
		with open(name, 'w') as f:
			for p in dictionaries:
				d = self.get_rosetta_energy_scores(p, models)
				first_key = d.keys()[0]
				second_key = d[first_key].keys()[0]
				header = d[first_key][second_key]["scores"][0].keys()
				header = ",".join(header)
				f.write("{0},{1},{2},{3},{4},{5},{6}\n".format(
					"model", "file_type", "analysis_type", "chain", "residue", "residue_number", header))
				for i in sorted(d):
					for j in sorted(d[i]):
						scores = [str(s)
								  for s in d[i][j]["scores"][0].values()]
						model = i.split(".")[0]
						chain = d[i][j]["chain"]
						residue_num = d[i][j]["residue_num"]
						residue = d[i][j]["residue"]
						string_of_scores = ",".join(scores)
						analysis_type = p
						super_string = "{0},{1},{2},{3},{4},{5},{6}\n".format(
							model, file, analysis_type, chain, residue, residue_num, string_of_scores)
						f.write(super_string)

	def output_sequence_logos(self, output, native, sequence="", prefix="", format="", units="", stacks_per_line="",
								title="", x_axis_label="", y_axis_height="", y_label="", errorbars="", fine_print="",
								color_scheme="", debug="",executable=""):
		
		seqlogos = SequenceLogos(output, native=native, designed=self.seq_analysis.designed_dict,
									sequence=sequence, output_prefix=prefix)
		
		command = seqlogos.initialize_options(format=format, units=units, stacks_per_line=stacks_per_line, title=title,
												y_axis_height=y_axis_height, y_label=y_label, errorbars=errorbars,
												fine_print=fine_print, color_scheme=color_scheme,path=executable,x_axis_label=x_axis_label)
		if debug:
			seqlogos.run_weblogo(command, debug=debug)
		else:
			seqlogos.run_weblogo(command)


class RosettaScores():

	def __init__(self, dict_of_interest, pdb):
		self.seq_analysis = dict_of_interest
		self.model = pdb
		self.score_table = scoretable.ScoreTable(pdb)

	def get_score(self):
		self.scores = {}
		for i in sorted(self.seq_analysis):
			res_score = self.score_table.get_score(chain=i[0], pdbres=i[1])
			self.scores[str(i[0]) + "_" + str(i[1])] = {"chain": i[0], "residue_num": int(
				i[1]), "scores": res_score, "residue": three_to_one(res_score[1])}
		return self.scores


class SequenceLogos():

	def __init__(self, pdbs, chain="", native=None, designed={}, sequence=True, output_prefix="Sequence_logo"):
		self.native = native
		self.pdbs = pdbs
		self.designed = designed.keys()
		self.chain = chain
		self.sequence = sequence
		self.output = output_prefix
		if self.native:
			if sequence == "sequence":
				self.annotation = self.get_native_sequence_and_numbering(
					self.native, self.designed)
			else:
				self.annotation = self.get_native_sequence_and_numbering(
					self.native, self.designed, sequence=self.sequence)

		if designed:
			for i in designed:
				self.designed.append(i)
			self.generate_fasta_from_pdbs_by_designed(self.pdbs, self.designed)
		else:
			self.generate_fasta_from_pdbs_by_chain(self.pdbs, self.chain)

	def generate_fasta_from_pdbs_by_chain(self, pdbs, chain_of_interest):
		'''takes in biopython pdb objects and outputs a fasta file, by chain'''
		with open(self.output + ".fasta", 'w') as f:
			for pdb in pdbs:
				for chain in pdb:
					if chain_of_interest == chain.id:
						this_chain = chain.id
						residues = []
						for residue in chain:
							residues.append(three_to_one(residue.resname))
						fasta_string = ">{0}_chain_{1}\n{2}\n".format(
							pdb.id, this_chain, "".join(residues))
						f.write(fasta_string)

	def generate_fasta_from_pdbs_by_designed(self, pdbs, designed):
		'''takes in biopython pdbs and outputs only the residues designed as a fasta'''
		with open(self.output + ".fasta", 'w') as f:
			for pdb in pdbs:
				residues = []
				for chain in pdb.get_chains():
					this_chain = chain.id
					for residue in chain:
						if (this_chain, residue.id[1]) in designed:
							residues.append(three_to_one(residue.resname))
				fasta_string = ">{0}_designed_residues\n{1}\n".format(
					pdb.id, "".join(residues))
				f.write(fasta_string)

	def get_native_sequence_and_numbering(self, native, design, sequence="sequence"):
		residues = []
		numbering = []
		for chain in native.get_chains():
			this_chain = chain.id
			for residue in chain:
				if (this_chain, residue.id[1]) in design:
					residues.append(three_to_one(residue.resname))
					numbering.append(this_chain + str(residue.id[1]))
		if sequence == "sequence":
			return residues
		else:
			return numbering

	def initialize_options(
		self, format="eps", sequence_type="protein", units="bits", first_index_number="1", lower_bound=None, upper_bound=None, size="large", stacks_per_line="40", title="Sequence Logo", x_axis_label=None, y_axis_height="4.32", y_label="",
		errorbars="No", fine_print="@Jordan Willis", color_scheme="classic", scale_width="Yes", composition="equiprobable",path=""):

		annotate = ",".join(self.annotation)
		# could of used Kwargs for all this s*it but just wrote it out
		self.input_file = ["-f", self.output + ".fasta"]
		self.seq_output = ["-o", self.output + "seq_log." + format]
		self.format = ["-F", format]
		self.sequence_style = ["-A", sequence_type]
		self.units = ["-U", units]
		self.first_index = ["-i", first_index_number]
		self.lower_bound = ["-l", lower_bound]
		self.upper_bound = ["-u", upper_bound]
		self.size = ["-s", size]
		self.stacks = ["-n", stacks_per_line]
		self.title = ["-t", title]
		self.annotate = ["--annotate", annotate]
		self.x_axis = ["-x", x_axis_label]
		self.y_axis_height = ["-S", y_axis_height]
		self.y_label = ["-y", y_label]
		self.errorbars = ["--errorbars", errorbars]
		self.color = ["-c", color_scheme]
		self.fine_print = ["--fineprint", fine_print]
		self.scale_width = ["--scale-width", scale_width]
		self.composition = ["--composition", composition]
		self.combined_arguments = self.input_file + self.seq_output + self.format + self.sequence_style + self.units + self.first_index + self.size + self.title + self.annotate +\
			self.y_axis_height + self.errorbars + self.color + \
			self.fine_print + self.scale_width + self.composition + self.stacks + self.y_label

		self.path = [path]
		if lower_bound:
			self.combined_arguments += self.lower_bound
		if upper_bound:
			self.combined_arguments += self.higher_bound
		if x_axis_label:
			self.combined_arguments += self.x_axis

		return self.combined_arguments

	def run_weblogo(self, command, debug=False):
		import subprocess as sub
		cmdline = self.path + command
		if debug:
			print " ".join(cmdline)
		sub.call(cmdline)
