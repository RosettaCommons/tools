import textwrap
import argparse
import sequence_analysis
import os.path

def args():
	# call on all our file type parsers in the sequence_anlysis_method
	type_parsers = sequence_analysis.DesignParsers()
	"""A customized argument parser that does a LOT of error checking"""
	parser = argparse.ArgumentParser(prog="design_analysis.py", formatter_class=argparse.RawTextHelpFormatter, description=textwrap.dedent('''\
									Design Analysis
					_____________________________________________________________________________________________\n
					This script is intented to encompass the entire functionality of design analysis.
					Everything you could want to do with design is called upon in this script.
					The most basic functionality is to pass a list of pdbs or get a position matrix of occurences count
					of just one line.

					The functionality extends from there by giving bitscores, changes in energy,
					giving position specific scoring matrices of your design,giving a customizable sequence logo. 
					This is a combination of many scripts and classes.\n
					\t\t\tJordan Willis'''))

	neces = parser.add_argument_group(
		title='Necessary', description="PDB files have to be included")
	neces.add_argument(
		"targets", metavar="*.pdb", nargs=argparse.REMAINDER,
		help="The PDB files to be analyzed")
	
	recommended = parser.add_argument_group(
		title='Recommended Options', description="Will give you a more complete analysis based on a res file, and a native pdb to compare it to.")
	recommended.add_argument("--native_pdb", "-p", dest="n_pdb",
						  help="The native pdb file to compare against", default=None)
	recommended.add_argument(
		"--corr", "-c", metavar="corr.corr", dest="corr", help="Get the results defined only in the corr file",
		type=type_parsers.parse_corr, default=None)
	recommended.add_argument(
		"--res", "-r", metavar="resfile.resfile", dest="rez", help="Get the results defined only in the residue file",
		type=type_parsers.parse_res, default=None)

	change_parser = parser.add_argument_group(title="Output Options: Please read carefully", description="These arguments change file name, which file is printed, which is output to a dictionary,and give verbose printing")
	change_parser.add_argument("--verbose", "-v", dest="verbose",
							   help="everything printed to a file will also be shown on the screen", action="store_true", default=False)
	change_parser.add_argument(
		"--prefix", "-P", dest="prefix", help="The prefix for what all the output files will be", default="output_")
	change_parser.add_argument("--score_files", "-s", dest="o_files", nargs="+", help="What do you want output to a file? Can list as a space seperated (eg -s n d nd):\n\
								a - full analysis dict\n\
								d - give analysis of just designed residues\n\
								n - just the native residues scores are shown\n\
								nd - just the native residues of the residues you designed\n\
								\nDefaults to full analysis dictionary", default=["a"])
	change_parser.add_argument(
		"-b", dest="bit", help="Should the output be in bit score? Defaults to occurences instead of bitscore", default=False,
							   action="store_true")
	change_parser.add_argument(
		"-S", dest="get_stats", help="If you specify a native file and a design file, it will give you an output of the stats of the design", 
								default=False, action="store_true")

	rosetta_parser = parser.add_argument_group(
		title="Rosetta Energy Analysis", description="Options for outputing options about energy scores,the dictionaries analyzed depend on what you asked for using the \
														-s output options flag")
	rosetta_parser.add_argument(
		"--rosetta", "-t", dest="term", default=False, action="store_true",
								help="This option will output a .csv file of the model,chain,residue,residue number,and rosetta scores.")
	
	bitscore = parser.add_argument_group(
		title="Bit Score Options'", description="options for bitscore metric for each designed residue")
	bitscore.add_argument(
		"-n", dest="normal", help="do you want the bit scores to be normalized by the shannon entropy",
		action="store_true", default=False)


	seqlogo_parser = parser.add_argument_group(
		title="Sequence Logos Options", description="These options handle the sequence logos that can be output from the design analysis script, and uses the api of weblogo to do so.")
	seqlogo_parser.add_argument("--seq", "-l", dest="logo", default=False, action="store_true",help="Turn on Sequence Logos for all the dictionaries you supply given in an .eps file")
	seqlogo_parser.add_argument("--path","-lp", dest="logo_path", default="/blue/meilerlab/apps/scripts_repository/sequence_analysis/weblogo-3.3/weblogo",help="What is the path to weblogo software? Defaults to meilerlab enviroment")
	seqlogo_parser.add_argument("--format", choices=['eps','jpg','png','png_print','pdf','jpeg','svg','logodata'],dest="s_format", default="eps",help="What format do you want the sequence logo in?")
	seqlogo_parser.add_argument("--units", choices=['bits','nats','kt','kJ/mol','kcal/mol','probability'],dest="s_units", default="bits",help="What do you want the units of the sequence logo to be in? Defaults to bits.")
	seqlogo_parser.add_argument("--stacks", dest="s_stacks", default="40",help="How many sequences per line in the logo, default=after forty letters it will go to a new line.")
	seqlogo_parser.add_argument("--title", dest="s_title", default="Sequence Logo",help="The title of your sequence logo")
	seqlogo_parser.add_argument("--x_label", dest="s_x_axis", default="",help="What do you want the x axis titled?")
	seqlogo_parser.add_argument("--y_axis_height", dest="s_y_height", default="4.32",help="How high do you want the Y axis, currently 4.32 which is the maximum acheivable score in a unbiased design")
	seqlogo_parser.add_argument("--y_label", dest="s_y_label", default="bits", help="Title of Y-Axis")
	seqlogo_parser.add_argument("--error_bars", dest="s_y_error", default="No",help="Do you want error bars turned on, YES/NO?")
	seqlogo_parser.add_argument("--fine_print", dest="s_fine",default="Jordan Willis", help="Fine Print")
	seqlogo_parser.add_argument("--color_scheme", choices=['auto','chemistry','charge','classic','hydrophobicity','monochrome'],dest="s_scheme", default="chemistry",help="The color scheme of the sequence logo. Defaults to Classic")
	seqlogo_parser.add_argument("--labels", dest="s_labels", choices=['sequece','numbers'], default="sequence", help="The x-axis labels can either take on the native residues sequence given with a native pdb file or the numbering of the pdb residue.")
	seqlogo_parser.add_argument("--debug", dest="s_debug", default=False, action="store_true",help="Get the full command line of what was put into weblogo")
	arg = parser.parse_args()

	target_type = ""

	# sometimes people will put the positional argumeent list right behind the
	# targets list, this is so we can differenciate between the two
	option_keys = ['a', 'd', 'n', 'nd']
	temp_tables = []
	for i in arg.o_files:
		if i in option_keys:
			temp_tables.append(i)
		else:
			arg.targets.append(i)

	arg.o_files = temp_tables

	native_pdb = arg.n_pdb
	target_list = arg.targets

	# check to see that their is one positional argument
	#@TODO, make these lamda functions in the argument parser
	if arg.n_pdb:
		try:
			arg.n_pdb = type_parsers.parse_pdbs(arg.n_pdb)
		except:
			parser.exit(
				status=1, message="\nnative pdb is not a pdb file or does not exist\n")

	if arg.targets == []:
		parser.print_help()
		parser.exit(
			status=1, message="\nNeed pdbs to analyze!!!!\ndesign_analysis.py *.pdbs\n")

	if ("d" in arg.o_files or "nd" in arg.o_files) and (arg.corr is None and arg.rez is None):
		parser.print_help()
		parser.exit(
			status=1, message="\n You have asked for a design dictionaires yet have not specified a .corr or .resfile to figure out what you designed\n")

	if ("n" in arg.o_files or "nd" in arg.o_files) and arg.n_pdb is None:
		parser.print_help()
		parser.exit(
			status=1, message="\n You have asked for a native dictionaires yet have not specified a native pdb to figure out what residues are native\n")
	try:
		# checking by iterating through the generator it returns and sees
		# if anything is there
		arg.targets = type_parsers.parse_pdbs(arg.targets)
		arg.targets[0].get_atoms().next()
		target_type = "pdb"
	except StopIteration:
		parser.print_help()
		raise TypeError(
			"The positional argument is not PDB")

	if arg.logo and arg.s_labels == "sequence" and arg.n_pdb is None:
		parser.exit(
			status=1, message="\nYou have asked for the native sequences to be on the x axis in the sequence logo but have not specified a native pdb\n")
	
	if arg.get_stats and arg.bit is False:
		parser.exit("You have asked for stats but did not ask for bitscore - that is currently not available.")

	if arg.logo and (not arg.rez and not arg.corr):
		 parser.exit("for seq logos, you must specifiy a res or corr file...can't possibly generate a seq logo of the entire PDB!")

	if arg.logo and ('n' not in arg.o_files and 'nd' not in arg.o_files):
		 parser.exit("You have asked for sequence logos but not specified design or native design output (eg. -s n or -s nd)")

	if arg.get_stats and 'nd' not in arg.o_files:
		parser.exit("You have asked for stats but didn't specify \"nd\" for the output type, which is what the stats are on")

	if arg.normal and not arg.bit:
		parser.exit("You have asked to normalize by shannon entropy but haven't asked for bit score, pass -b")
	
	if arg.logo:
		if not os.path.isfile(arg.logo_path):
			parser.exit("The sequence logo execution file {0} you have specified does not exist".format(arg.logo_path))
	return arg, target_type, native_pdb, target_list, native_pdb