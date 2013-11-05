#!/usr/bin/env python2.7

import sys
if sys.version_info < (2, 7):
	raise Exception("You must use python2.7 to run this")


import sequence_analysis as sequence_analysis
import analysis_tools as at
from argument_parsing import args

def main():

	# get back all the arguments that have been parsed
	arg, target_type, native_pdb, target_list, native_pdb = args()
	# call on the class DesignScoreTable That does the dirty work. Initially I
	# was going to analyze everything, but I think I want to call on
	# dictionaries as I get them
	design_score_table = sequence_analysis.DesignScoreTable(native_pdb=arg.n_pdb,
											 t_type=target_type,
											 targets=arg.targets,
											 corr_file=arg.corr,
											 rez_file=arg.rez,
											 out_files=arg.o_files)

	#number of pdb files
	pdb_file_number = len(arg.targets)
	analyzer_t = at.AnalysisTools(
		design_score_table,pdb_file_number)  # I on the other hand would like to anlysis tools, a clean interface that does not show the backend.

	# print verbose dictionaries defined in o files
	if arg.verbose:
		#print out dictionary files verbosely, should be decorator
		analyzer_t.verbose(arg.o_files)

	# output dictionaries to files
	
	analyzer_t.output(
		prefix=arg.prefix, list_o=arg.o_files, bit=arg.bit,normal=arg.normal)

	if arg.get_stats:
		analyzer_t.get_stats(filename=arg.prefix + "_stats.txt",bit=arg.bit)

	# get rosetta energy scores
	if arg.term:
		analyzer_t.output_rosetta_dict(arg.o_files, arg.prefix, target_list)

	#finally, run sequence logo
	if arg.logo:
		 analyzer_t.output_sequence_logos(
			arg.targets, arg.n_pdb, format=arg.s_format, units=arg.s_units, stacks_per_line=arg.s_stacks, title=arg.s_title, x_axis_label=arg.s_x_axis,
				y_axis_height=arg.s_y_height, y_label=arg.s_y_label, errorbars=arg.s_y_error, fine_print=arg.s_fine, color_scheme=arg.s_scheme, 
				sequence=arg.s_labels, prefix=arg.prefix, debug=arg.s_debug,executable=arg.logo_path)
