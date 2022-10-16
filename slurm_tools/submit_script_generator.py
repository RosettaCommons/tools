#!/usr/bin/python
"""
Make folders, flags, sbatch files for MARCC submission for LK peptide surface 
docks. Vary weight parameters, surface atom properties, hydrophobic/hydrophilic
surface.
Created by Joseph Lubin
"""

import argparse
import os 
import numpy as np
import math

# Collecting inputs
parser = argparse.ArgumentParser()
subparser = parser.add_subparsers(help = "How many parameters to vary", 		\
	dest = 'generator_type')

# Options
parser.add_argument("-f", "--subfolder", help = 								\
	"What folder below the current working directory do you want to use?")
parser.add_argument("-s", "--score", help = 									\
	"What score function do you want to use, if not talaris2013? (beta_nov15)")
parser.add_argument("-c", "--capped", action = "store_true", help = 			\
	"Should the peptides be N-acetylated and C-amidated, as Collier & Latour?")
parser.add_argument("-n", "--n_cap", action = "store_true", help = 				\
	"Should the peptides be N-acetylated, as Weidner & Castner? (14mers only)")
parser.add_argument("-d", "--decoys", type = int, default = 10000, help = 		\
	"How many decoys do you want? (Default = 10000)")
parser.add_argument("-r", "--range", type = float, default = 5, help = 			\
	"Desired range cutoff for LJ & LK potentials? (Default = 5)")
parser.add_argument("-kcd", type = float, default = 0.52, help = 				\
	"Altered hydrophobicity LK_DGFREE of Lysine delta carbon? <0=hydrophilic")
parser.add_argument("-kce", type = float, default = 0.52, help = 				\
	"Altered hydrophobicity LK_DGFREE of Lysine epsilon carbon? <0=hydrophilic")
parser.add_argument("-e", "--exclusions", type = str, action = 'append', 		\
	choices = "abio74", help = "What test conditions do you want to exclude? " +\
	"Can be used multiple times to exclude multiple categories. Sequence: " + 	\
	"a--alpha, b--beta // Surface: i--" + "hydrophilic, o--hydrophobic // " +	\
	"Size: 7--7mer, 4--14mer")
parser.add_argument("-x", "--execute", action = "store_true", help =  			\
	"Do you wish to submit all generated scripts?")

	#Baseline case
base = subparser.add_parser('b', help = "baseline")

	#Single parameter case
single = subparser.add_parser('s', help = "Vary one parameter")	
single.add_argument("parameter_1", type = str, help =  							\
	"What parameter do you mean to vary?")
single.add_argument("begin_1", type = float, help = 							\
	"Minimum of desired parameter range for testing")
single.add_argument("end_1", type = float, help = 								\
	"Maximum of desired parameter range for testing")
single.add_argument("step_1", type = float, help = "Desired increments in range")

	#Double parameter case
double = subparser.add_parser('d', help = "Vary two parameters")
p1 = double.add_argument_group("first parameter")
p1.add_argument("parameter_1", type = str, help =  								\
	"What is the first parameter you mean to vary?")
p1.add_argument("begin_1", type = float, help = 								\
	"Minimum of first parameter range for testing")
p1.add_argument("end_1", type = float, help = 									\
	"Maximum of first parameter range for testing")
p1.add_argument("step_1", type = float, help = 									\
	"Desired increments in first parameter range")

p2 = double.add_argument_group("second parameter")
p2.add_argument("parameter_2", type = str, help = 								\
	"What is the second parameter you mean to vary?")
p2.add_argument("begin_2", type = float, help = 								\
	"Minimum of second parameter range for testing")
p2.add_argument("end_2", type = float, help = 									\
	"Maximum of second parameter range for testing")
p2.add_argument("step_2", type = float, help = 									\
	"Desired increments in second parameter range")

args = parser.parse_args()

#################################################################################
#Defining parameter lists and functions
def get_parameter_type(parameter):
	"""
	Given a parameter name, this function checks the Rosetta default weight 
	parameters and atom properties, and matches the type of the parameter
	"""
	weight_params = ['fa_atr', 'fa_rep', 'fa_sol', 'fa_intra_rep', 'fa_elec', 	\
					'pro_close', 'hbond_sr_bb','hbond_lr_bb', 'hbond_bb_sc', 	\
					'hbond_sc', 'dslf_fa13', 'rama','omega', 'fa_dun', 			\
					'p_aa_pp', 'ref']
	atom_properties = ['LJ_RADIUS', 'LJ_WDEPTH', 'LK_DGFREE', 'LK_LAMBDA', 		\
					'LK_VOLUME']
	parameter_types_list = ["Weight", "Atom Property", "Lysine CD & CE", "ERROR"]

	if parameter.lower() in weight_params:
		par_type = "weight"
		print_index = 0 

	elif parameter.upper() in atom_properties:
		par_type = "atom_property"
		print_index = 1
		print "\n\nParameter type: "

	elif parameter.lower() == 'kcde':
		par_type = "kcde"
		print_index = 2

	else: 
		par_type = "ERROR"
		print_index = 3

	print parameter + ":\t" + parameter_types_list[print_index]

	return par_type

def trial_conditions_chooser(exclusions):
	"""
	This function takes an input with letter or number codes indicating a data 
	set to exclude from testing.  Options are alpha or beta peptide,hydrophobic 
	or hydrophilic surface, and 7mer or 14mer peptide.
	"""
	condition_files =  ["hydrophilic_14mer_alpha", "hydrophobic_14mer_alpha", 	\
						"hydrophilic_14mer_beta", "hydrophobic_14mer_beta", 	\
						"hydrophilic_7mer_alpha", "hydrophobic_7mer_alpha", 	\
						"hydrophilic_7mer_beta", "hydrophobic_7mer_beta"]

	test_cond = []

	seq_exclude = [0,1]
	surf_exclude = [0,1]
	size_exclude = [0,1]

	if "14" in exclusions:
		del size_exclude[0]
	elif "7" in exclusions or args.n_cap:
		del size_exclude[1]

	if "a" in exclusions:
		del seq_exclude[0]
	elif "b" in exclusions:
		del seq_exclude[1]

	if "i" in exclusions:
		del surf_exclude[0]
	elif "o" in exclusions:
		del surf_exclude[1]

	for i in size_exclude:
		for j in seq_exclude:
			for k in surf_exclude:
				file_index = 4*(i)+2*(j)+(k)
				test_cond.append(condition_files[file_index])

	test_cond_display = str(test_cond)
	text_replacements = {"[":"", "]":"", "'": "", "_": " "}
	for i, j in text_replacements.iteritems():
		test_cond_display = test_cond_display.replace(i, j)

	print "Test conditions:"
	print test_cond_display
	return test_cond

def freerange(start,stop,step):
	"""
	Generates a range that doesn't need to count by 1
	"""
	i = start
	while round(i, 3) <= stop:
		yield round(i, 3)
		i += step

def range_gen(start, stop, steps):
	"""
	Uses the freerange function to generate a list of parameter variant values in 
	both numerical and text form
	"""
	param_set = []
	variable_step = (stop - start)/steps

	for i in freerange(start, stop, variable_step):
		param_set.append(i)

	return param_set

def determine_score_file(parameter_value_1, parameter_1, parameter_type_1, 		\
						parameter_value_2, parameter_2, parameter_type_2):
	""" Determines appropriate score weights file to flag"""
	param_list = {'fa_atr': 0.8, 'fa_rep': 0.44, 'fa_sol': 0.75, 'fa_intra_rep':\
	0.004, 'fa_elec': 0.7, 'pro_close': 1, 'hbond_sr_bb': 1.17, 'hbond_lr_bb': 	\
	1.17, 'hbond_bb_sc': 1.17,'hbond_sc': 1.1,'dslf_fa13': 1.0,	'rama': 0.2, 	\
	'omega': 0.5, 'fa_dun': 0.56, 'p_aa_pp': 0.32, 'ref': 1}

	# Pointing to appropriate scores file
	if args.score:
		scorefile = args.score
	else:
		if [parameter_type_1, parameter_type_2].count("weight") == 0:
			scorefile = 'talaris2013'
		elif [parameter_type_1, parameter_type_2].count("weight") == 1:
			if parameter_type_1 == "weight":
				if parameter_value_1 == param_list[parameter_1.lower()]:
					scorefile = 'talaris2013'
				else:
					scorefile = 'talaris_' + parameter_1.lower() + '_' + 		\
						str(parameter_value_1)
			else:
				if parameter_value_2 == param_list[parameter_2.lower()]:
					scorefile = 'talaris2013'
				else:
					scorefile = 'talaris_' + parameter_2.lower() + '_' + 		\
					str(parameter_value_2)
		else:
			if parameter_value_1 == param_list[parameter_1.lower()] and 		\
						parameter_value_2 == param_list[parameter_2.lower()]:
				scorefile = 'talaris2013'
			else:
				scorefile = 'talaris_'+ parameter_1.lower() + '_' + 			\
					str(parameter_value_1) + '_' + parameter_2.lower() + '_' + 	\
					str(parameter_value_2)

	return scorefile

def kcde_value_flags(parameter_value):
	""" 
	Sets KCD value based on KCE value. Changes by 14.4% of KCE, so KCD is 1 when 
	KCE is 10.
	"""
	kcde_lines = []
	atom_adjust = '-chemical:set_atom_properties fa_standard:'

	if parameter_value == 'manual':
		kcd_value = args.kcd
		kce_value = args.kce
	else: 
		kce_value = parameter_value
		kcd_value = 0.52 - round(0.144 * (0.52 - kce_value),2)

	kcde_lines.append(atom_adjust + 'KCD:LK_DGFREE:'+ str(kcd_value) + '\n')
	kcde_lines.append(atom_adjust + 'KCE:LK_DGFREE:'+ str(kce_value) + '\n')

	return kcde_lines

def flags_maker(path, destination, rosetta_folder, source_folder, dataset, 		\
				num_decoys, detection_range, parameter_value_1, parameter_1, 	\
				parameter_type_1, parameter_value_2, parameter_2, 				\
				parameter_type_2):
	"""
	Writes flags file for LK peptide SAM surface docking simulation
	"""
	loc_ident = os.path.join(path, destination, destination)
	flagname = destination + ".flags"
	scorefile = determine_score_file(parameter_value_1, parameter_1, 			\
		parameter_type_1, parameter_value_2, parameter_2, parameter_type_2)
	
	#addressing kcde values
	kcde_lines = []
	if parameter_1 == 'kcde':
		kcde_lines = kcde_value_flags(parameter_value_1)
	elif parameter_2 == 'kcde':
		kcde_lines = kcde_value_flags(parameter_value_2)
	elif args.kcd != 0.52 or args.kce != 0.52:
		kcde_lines = kcde_value_flags('manual')

	# Writing flags file
	with open(flagname,'w') as fl:
		fl.write('-database ' + rosetta_folder + 'database\n')
		fl.write('-include_surfaces\n')
		if args.capped:
			fl.write('-s ' + source_folder + dataset + '_capped.pdb\n')
		elif args.n_cap:
			fl.write('-s ' + source_folder + dataset + '_n_capped.pdb\n')
		else:
			fl.write('-s ' + source_folder + dataset + '.pdb\n')
		fl.write('-in:file:surface_vectors ' + source_folder + 'SAM.surf\n')
		fl.write('-mute core\n')
		fl.write('-mute protocols.moves.RigidBodyMover\n')
		fl.write('-nstruct ' + num_decoys + '\n')
		if args.range != 5:
			fl.write('-score::fa_max_dis ' +  str(detection_range) + '\n')
		fl.write('-score:weights ' + scorefile + '\n')
		if parameter_type_1 == "atom_property":
			fl.write('-chemical:set_atom_properties fa_standard:CH3S:' + 		\
				parameter_1.upper() + ':' + str(parameter_value_1) + '\n')
		if parameter_type_2 == "atom_property":
			fl.write('-chemical:set_atom_properties fa_standard:CH3S:' + 		\
				parameter_2.upper() + ':' + str(parameter_value_2) + '\n')
		for line in kcde_lines:
			fl.write(line)
		fl.write('-out:pdb_gz\n')
		fl.write('-out:path:pdb '+ loc_ident + '_output\n')
		fl.write('-out:path:score '+ loc_ident + '_output\n')
		fl.write('-mpi_tracer_to_file ' + loc_ident + '_err/tracer.out\n')

	os.rename(flagname, os.path.join(path, destination, flagname))
	print "\t" + flagname
	return flagname

def sbatch_maker(flags_file, path, destination, rosetta_folder, condition_name, \
				wall_time, parameter_name_1, set_member_1, parameter_name_2, 	\
				set_member_2):
	"""
	Writes MARCC sbatch file for LK peptide SAM surface docking simulation
	"""
	loc_ident = os.path.join(path, destination, destination)
	subname = destination+".sbatch"

	# Dictionaries to abbreviate conditions and parameters
	conditions_codes = {"hydrophobic_14mer_alpha": "O4A", 						\
						"hydrophilic_14mer_alpha": "I4A", 						\
						"hydrophobic_14mer_beta": "O4B", 						\
						"hydrophilic_14mer_beta": "I4B", 						\
						"hydrophobic_7mer_alpha": "O7A", 						\
						"hydrophilic_7mer_alpha": "I7A", 						\
						"hydrophobic_7mer_beta": "O7B", 						\
						"hydrophilic_7mer_beta": "I7B"}

	parameter_codes = {'fa_atr': 'ATR', 'fa_rep': 'REP', 'fa_sol': 'SOL', 		\
					'fa_intra_rep': 'IRP', 'fa_elec': 'ELC', 'pro_close': 'PCL',\
					'hbond_sr_bb': 'HSB', 'hbond_lr_bb': 'HLB', 				\
					'hbond_bb_sc': 'HBS', 'hbond_sc': 'HSC', 'dslf_fa13': 'DSF',\
					'rama': 'RAM', 'omega': 'OMG', 'fa_dun': 'DUN', 			\
					'p_aa_pp': 'PAP', 'ref': 'REF', 'LJ_RADIUS': 'RAD', 		\
					'LJ_WDEPTH': 'WDP', 'LK_DGFREE': 'DGF', 'LK_LAMBDA': 'LAM', \
					'LK_VOLUME': 'VOL', 'kcde': 'KC','baseline': 'BAS', '': ''}

	# Getting parameter type
	try:
		par_code_1 = parameter_codes[parameter_name_1]
	except:
		par_code_1 = parameter_codes[parameter_name_1.upper()]

	try:
		par_code_2 = parameter_codes[parameter_name_2]
	except:
		par_code_2 = parameter_codes[parameter_name_2.upper()]

	# Writing SBATCH files
	with open(subname,'w') as su:
		su.write('#!/bin/bash -l\n')
		su.write('\n')
		su.write('#SBATCH --job-name=' + conditions_codes[condition_name] + 	\
			set_member_1 + par_code_1 + set_member_2 + par_code_2 + '\n')
		su.write('#SBATCH --partition=parallel\n')
		su.write('#SBATCH --time=' + wall_time + ':0' + "\n")
		su.write('#SBATCH --nodes=5\n')
		su.write('#SBATCH --ntasks-per-node=24\n')
		su.write('#SBATCH --mem=120GB\n')
		su.write('#SBATCH --output ' + loc_ident + '_err/report.%j.out\n')
		su.write('#SBATCH --error ' + loc_ident + '_err/report.%j.err\n')
		su.write('#SBATCH --mail-user=jlubin3@jhu.edu\n')
		su.write('#SBATCH --mail-type=ALL\n')
		su.write('\n')
		su.write('module unload openmpi\n')
		su.write('module load intel-mpi\n')
		su.write('\n')
		su.write('ROSETTABIN=' + rosetta_folder + 'source/bin\n')
		su.write('ROSETTAEXE=surface_docking\n')
		su.write('COMPILER=mpi.linuxgccrelease\n')
		su.write('EXE=$ROSETTABIN/$ROSETTAEXE.$COMPILER\n')
		su.write('echo Starting MPI job running $EXE\n')
		su.write('\n')
		su.write('date\n')
		su.write('time mpirun $EXE @' + path + '/' + destination + '/' + 		\
			flags_file + '\n')
		su.write('date\n')
	os.rename(subname, os.path.join(path, destination, subname))
	print "\t" + subname
	return subname, loc_ident

def par_value_stringer(value):
	""" 
	When given a value, this function will convert it into a string, with a 
	consistent length so folders sort nicely.
	"""
	string = str(value)
	if '.' not in string:
		string += '.0'
	while len(string) < 5:
		string += '0'

	return string

def marcc_file_maker(path, test_condition, decoy_count, detection_range, 		\
					set_item_1, set_item_2, parameter_1, parameter_value_1, 	\
					parameter_type_1, parameter_2, parameter_value_2, 			\
					parameter_type_2):
	"""
	Uses flags_maker and sbatch_maker and generates flags, sbatch scripts,
	and folders for LK peptide surface docking simulations run on MARCC.
	"""
	# making folder
	if parameter_type_1 == "base":
		newpath = test_condition + "_baseline"
	elif parameter_type_2 == "":
		newpath = test_condition + "_" + parameter_1 + "_" + 					\
			par_value_stringer(parameter_value_1)
	else:
		newpath = test_condition + "_" + parameter_1 + "_" + 					\
			par_value_stringer(parameter_value_1) + "_" + parameter_2 + "_" + 	\
			par_value_stringer(parameter_value_2)

	os.makedirs(os.path.join(path, newpath))
	print "Folder: " + newpath

	ros_folder = "/home-2/jlubin3@jhu.edu/work/jgray21/jhlubin/Rosetta/main/"
	work_folder = "/home-2/jlubin3@jhu.edu/scratch/" + 							\
		"samonolayer_parameter_refine/seed_files/"

	# making flags
	flagname = flags_maker(path, newpath, ros_folder, work_folder, 				\
						test_condition, decoy_count, detection_range, 			\
						parameter_value_1, parameter_1, parameter_type_1, 		\
						parameter_value_2, parameter_2, parameter_type_2)

	# making submission script
	subname, loc_ident = sbatch_maker(flagname, path, newpath, ros_folder, 		\
									test_condition, time, parameter_1, 			\
									set_item_1, parameter_2, set_item_2)

	subscript = (path + "/" + newpath + "/" + subname)

	#making output and error folders
	outname = loc_ident + "_output"
	os.makedirs(outname)
	print "\tFolder: " + newpath + "_output"
	ername= loc_ident + "_err"
	os.makedirs(ername)
	print "\tFolder: " + newpath+ "_err\n"

	return subscript, subname

#################################################################################
print "\n"

# Checking parameter type
if args.generator_type == "b":
	par_type_1 = "base"
	parameter_1 = "baseline"
	par_type_2 = ""
	parameter_2 = ""
	print "Parameter type: "
	print "Baseline"

elif args.generator_type == "s":
	print "Parameter type: "
	par_type_1 = get_parameter_type(args.parameter_1)
	parameter_1 = args.parameter_1
	par_type_2 = ""
	parameter_2 = ""

elif args.generator_type == "d":
	print "Parameter type: "
	par_type_1 = get_parameter_type(args.parameter_1)
	par_type_2 = get_parameter_type(args.parameter_2)	
	parameter_1 = args.parameter_1
	parameter_2 = args.parameter_2
	if par_type_2 == "ERROR":
		par_type_1 = par_type_2

print ""

#Response to incorrect parameter
if par_type_1 == "ERROR":
	print "\nPlease check the name of the parameter you want to vary."
	print "\tWeight parameters:"
	for i in range(len(weight_params)):
		print "\t\t" + weight_params[i]
	print "\tAtom Properties:"
	for i in range(len(atom_properties)):
		print "\t\t" + atom_properties[i]
	print "\tLysine CD & CE:\n\t\tkcde"
	exit()

# Data sets to generate
if args.exclusions:
	test_cond = trial_conditions_chooser(args.exclusions)
else:
	test_cond = trial_conditions_chooser([])

# Generating list of values from given min, max, and step count
if args.generator_type == "b":
	param_set_1 = [0]
	param_set_2 = [0]

elif args.generator_type == "s":
	param_set_1 = range_gen(args.begin_1, args.end_1, args.step_1)
	param_set_2 = [0]
	print "\nList of parameter values:"
	print param_set_1

elif args.generator_type == "d":
	param_set_1 = range_gen(args.begin_1, args.end_1, args.step_1)
	param_set_2 = range_gen(args.begin_2, args.end_2, args.step_2)
	print "\nList of parameter values:"
	print parameter_1 + ":\t" + str(param_set_1)
	print parameter_2 + ":\t" + str(param_set_2)

# Adjusting time for increased for decoy counts
decoys = str(args.decoys).replace(".0", "")
hours = math.ceil(args.decoys / 1000) * 1.5 + 0.5

if hours % 1 == 0:
	hours = str(hours).replace(".0", "")
	minutes = "0"
else:
	minutes = math.ceil((hours % 1) * 60)
	minutes = str(minutes).replace(".0", "")
	hours = math.floor(hours)
	hours = str(hours).replace(".0", "")
time = hours + ":" + minutes

print "\n" + decoys + " decoys"
print "Requesting time per simulation:\t" + time

trials = len(test_cond) * len(param_set_1) * len(param_set_2) * int(decoys)
print "Total decoys generated:\t" + str(trials)

# Setting destination
if args.subfolder:
	path = os.path.join(os.getcwd(), args.subfolder)
	if not os.path.isdir(path):
		os.makedirs(path)
		print "\nCreated destination folder:"
	else:
		print "Destination folder:"	
	print path		
else: 
	path = os.getcwd()
	print "\nDestination folder:"
	print path

# Creating files
subscripts = []
subnames = []
print "\n\nMaking the following:"
for i in range(len(test_cond)):
	for j in range(len(param_set_1)):
		for k in range(len(param_set_2)):
			subscript, subname = marcc_file_maker(path, test_cond[i], decoys, 	\
										args.range, str(j), str(k), 			\
										parameter_1, param_set_1[j], par_type_1,\
										parameter_2, param_set_2[k], par_type_2)

			subscripts.append(subscript)
			subnames.append(subname)

# Submitting the sbatch files
#Should only be done after the weight files are generated, if varying a weight
if args.execute:
	for i in range(len(subscripts)):
		bashline = "sbatch " + str(subscripts[i])
		os.system(bashline)
	print "\nAll jobs submitted"

# Printing list of submission scripts
#The conversion below makes it easier to execute the list of subscripts with a 
#loop bash command (for i in [block]; do sbatch $i; done) if the execute 
#option is not used
else:
	sub_string = ""
	for i in subscripts:
		sub_string = sub_string + i
		if i != len(subscripts):
			sub_string = sub_string + " "
	print "\nList of sbatch files:"
	print sub_string
	print "\n"
