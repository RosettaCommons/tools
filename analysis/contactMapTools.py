#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
# @brief    utility for performing several tasks with contactMaps
# @detail
# @author   Joerg Schaarschmidt (schaarscjoe@medizin.uni-leipzig.de)
# @version  1.0
# @date     07/25/2012

import sys
import argparse
import math
import re
from types import *

# Define Global Variables
float_digits = 3

# Main function - Argument parsing and call of subcommands ---------------------
def main():
	# Process commandline arguments - use -h for usage statement
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--aacode',dest='aacode', type=int, default=0, help="Output Residue names in three (3) or one (1) letter code")
	parser.add_argument('-o', '--outfile',dest='outfile', default='', help="Output filename (default: stdout)")
	subparsers = parser.add_subparsers(help="type »./contactMapTools.py <subcommand> -h« for help regarding the subcommands")

	# create the parser for the "combine" command
	parser_combine = subparsers.add_parser('combine',description="%(prog)s combines multiple contact map files (matrix format) by adding up the contact counts of the single files")
	parser_combine.add_argument('files', nargs='*', help="the contact maps to combine")
	parser_combine.add_argument('-l','--filelist',help="list of contact map files to combine in plain text format", type=argparse.FileType('r'))
	parser_combine.add_argument('-f', '--fraction', help="print the fraction of models with the contact instead of the model count ",action="store_true", default=False)
	parser_combine.set_defaults(func=combine)

	# create the parser for the "compare" command
	parser_compare = subparsers.add_parser('compare',description="%(prog)s compares the contact map of a structure or set of structures to the contact map of a single reference structure.")
	parser_compare.add_argument('native_file', help="reference map of a single structure to which the second map will be compared")
	parser_compare.add_argument('model_file', help="map of the structure(s) that will be compared to the reference map")
	parser_compare.add_argument('-s','--seqsep', type=int, default=0, help="minimal distance that two residues have to be seperated in the amino acid sequence to be included in the comparison [def: 0]")
	parser_compare.add_argument('-c','--cutoff', type=float, default=0.1, help="fraction or number of models the contact count has to exceed for a contact to be considered positive [def: 0.1]")
	parser_compare.set_defaults(func=compare)

	# create the parser for the "roc_curve_analysis" command
	parser_roc = subparsers.add_parser('roc',description="%(prog)s compares the contact map of a structure or set of structures to the contact map of a single reference structure and calculets the parameters for ROC analysis based on frequency of the contact in the set of models.")
	parser_roc.add_argument('native_file', help="reference map of a single structure to which the second map will be compared")
	parser_roc.add_argument('model_file', help="map of the structure(s) that will be compared to the reference map")
	parser_roc.add_argument('-s','--seqsep', type=int, default=0, help="minimal distance that two residues have to be seperated in the amino acid sequence to be included in the comparison [def: 0]")
	parser_roc.add_argument('-e','--exclude',help="Exclude biased contacts (present in all or none of the models", action="store_true", default=False)
	parser_roc.add_argument('-t','--statistics',help="output additional statistics for ROC-curve analysis", action="store_true", default=False)
	parser_roc.add_argument('-p','--pdbid', default="", help="optional PDB-ID or tag to add as additional column")
	parser_roc.set_defaults(func=roc_curve_analysis)


	# create the parser for the "renumber" command
	parser_renumber = subparsers.add_parser('renumber',description="Renumber the row and/or column names of a contact map (default matrix format) NOTE: input files will be modified !!!")
	parser_renumber.add_argument('input_files', nargs='+', help="Files to be renumbered")
	parser_renumber.add_argument('-p','--position',dest='pos',type=str,default='both',choices=['row', 'col', 'both'],help="Renumber only row names [row], column names [col], or both [both]")
	parser_renumber.add_argument('-s','--start', metavar='START(,COL_START)',required=True, help="New Start of Residue numbers in format <start> e.g. '15' or separate for rows and columns <rowstart>,<colstart> e.g '1,27'")
	parser_renumber.add_argument('-o', '--offset' , default='0', help="Start renumbering at entry i of row names and j of col names (format like -s option [i] or [i,j])")
	parser_renumber.set_defaults(func=renumber)

	parser_convert = subparsers.add_parser('convert',description="Convert Residue names in three (3) or one (1) letter code NOTE: input files will be modified !!!")
	parser_convert.add_argument('input_files', nargs='+', help="Files to be Converted")
	parser_convert.add_argument('-a', '--aacode',dest='aacode', type=int, default=1, help="Output Residue names in three (3) or one (1) letter code")
	parser_convert.set_defaults(func=convert)

	# parse the args and call whatever function was selected
	options = parser.parse_args()
	options.func(options)


# Shared functions (alphabetical order) ----------------------------------------

# Add contacts of the contact map to the combined contact map
def add_contact_map(file, matrix, models):
	nmatrix, nmodels = read_contact_map_file(file, True)
	for row in range(1,len(matrix)):
		for col in range(1,len(matrix[1])):
			# Increase the contact count of the current cell by the corresponing
			# contact count of the curent file
			matrix[row][col] += nmatrix[row][col]
	return (matrix, models+nmodels)


# Convert 3 to 1 letter aa code and vice versa
def convert_aa_names(matrix, aacode):
	if aacode != 0 :
		class cdict(dict):
			def __missing__(self, key):
				return key

		# dictionary for three to one-letter code conversion

		aa3to1=cdict(ALA='A', CYS='C', ASP='D', GLU='E', PHE='F', GLY='G',
			HIS='H', ILE='I', LYS='K', LEU='L', MET='M', ASN='N', PRO='P',
			GLN='Q', ARG='R', SER='S', THR='T', VAL='V', TRP='W', TYR='Y')
		# dictionary for one to three-letter code conversion
		aa1to3=cdict(A='ALA', C='CYS', D='ASP', E='GLU', F='PHE', G='GLY',
			H='HIS', I='ILE', K='LYS', L='LEU', M='MET', N='ASN', P='PRO',
			Q='GLN', R='ARG', S='SER', T='THR', V='VAL', W='TRP', Y='TYR')

		def convertName(residue_string):
			split_pos=re.search('[0-9]',residue_string).start(0)
			if split_pos==3 and aacode ==1 :
				return aa3to1[residue_string[:split_pos]]+residue_string[split_pos:]
			elif split_pos==1 and aacode ==3:
				return aa1to3[residue_string[:split_pos]]+residue_string[split_pos:]
			else:
				return residue_string

		for col in range (1, len(matrix[0])):
			matrix[0][col]= convertName(matrix[0][col])
		for row in range(1, len (matrix)):
			matrix[row][0]= convertName(matrix[row][0])

	return matrix


# Define function for parsing option strings
def parse_option_string(opt_string='1,2'):
	split_pos =opt_string.find(',')
	if split_pos==-1 :
		return [int(opt_string), int(opt_string)]
	else :
		return [int(opt_string[:split_pos]), int(opt_string[split_pos+1:])]


# Print ContactMap
def print_map(matrix, header, filename=''):
	fh = sys.stdout if filename =='' else open(filename, 'w')
	fh.write(header+'\n\n')
	for row in matrix:
		output_string= str(row[0])
		for col in row[1:]:
			if type(col) is FloatType:
				value_string = '{:.{}f}'.format(float(col), float_digits)
			else:
				value_string = str(col)
			output_string += '\t' + value_string
		fh.write(output_string+'\n')
	fh.flush()
	if  filename !='':
		fh.close()


# Read in contact map
def read_contact_map_file(filename , castInt = False):
	# Initialize contact matrix
	matrix=[]
	models=0
	comments = ''
	# Loop over lines of the current file
	for line in open(filename, 'r'):
		# Remove newline and trailing whitespace
		line=line.rstrip()
		# Extract the number of models from the comment lines
		if re.match("^#", line) :
			if re.match("^ Number of Models:", line):
				words = line.split()
				models += int (words[4])
			else:
				comments += line + '\n'
		# Split current row into columns and append to matrix
		elif line:
			matrix.append(line.split('\t'))
	if castInt:
		# Loop over all the cells of the matrix and cast the contact counts to
		# int values
		for row in range(1,len(matrix)):
			for col in range(1,len(matrix[1])):
				# cast contact count of current cell to int
				matrix[row][col] = int (matrix[row][col])
	return (matrix, models, comments)


# Convert the contact count to a percantage value
def convert_to_fraction(matrix, models):
	for row in range(1,len(matrix)):
		for col in range(1,len(matrix[1])):
			if type(col) is IntType :
				matrix[row][col] = float(matrix[row][col])/models
	return matrix


# Combine Subcommand -----------------------------------------------------------
# @brief    combines multiple contact map files
# @detail   combines multiple contact map files(only matrix format) by adding up
#           the contact counts and prints out a combined map to STDOUT
def combine(options):
	# Initialize variables
	first_file = ''
	# Retrieve filenme of first file to process
	if (options.filelist):
		first_file = options.filelist.readline().rstrip()
	else:
		first_file =options.files[0]

	# Read in first file
	matrix, models, comments = read_contact_map_file(first_file, True)

	# Loop over the remaining files and add them to the combined contact map
	filenames = options.filelist.readlines().rstrip() if options.filelist else options.files[1:]
	for file in filenames:
		matrix, models = add_contact_map(file, matrix, models)

	# Convert the contact count to a percantage value if -f option was used
	if (options.fraction):
		matrix=convert_to_fraction(matrix, models)

	matrix = convert_aa_names(matrix, options.aacode)

	# Print the combined map
	header='# Number of Models: {}'.format(models)
	print_map(matrix, header, options.outfile)


# Compare subcommand -----------------------------------------------------------
# @brief    compares two contact maps
# @detail   compares the contact map of a structure or set of structures to
#           the contact map of a single reference structure and depending on
#           the cutoff outputs a map highlighing the similarities and
#           differences of the two maps including a summary
#           optionally a roc curve for the cutoff value can be generated
def compare(options):

	# Initialize variables
	matrix, models, comments = read_contact_map_file(options.native_file, True)
	matrix_mod, models, comments = read_contact_map_file(options.model_file, True)
	nrow=0
	nignored=0
	ntruepos=0
	nfalsepos=0
	ntrueneg= 0
	nfalseneg= 0
	for row in range(1,len(matrix)):
		for col in range(1,len(matrix[1])):
			# Skip contacts that violate specified sequence separation
			if (math.fabs(col-row)<options.seqsep):
				nignored +=1
				matrix[row][col] = 0.0 # unconsidered residue ->0
			# Check if value of Contact Maps qualifies as positive value
			elif (float(matrix_mod[row][col])>=options.cutoff):
				# Compare positive entry to reference structure and assign value accordingly
				if (matrix[row][col] == 1):
					matrix[row][col] = 1.0 # true positive -> 1
					ntruepos += 1
				else:
					matrix[row][col] = +0.5 # false positive -> 0.5
					nfalsepos += 1

			else:
				# Compare nagative entry to reference structure  and assign value accordingly
				if (matrix[row][col] == 1):
					matrix[row][col] = -0.5 # false negative -> -0.5
					nfalseneg+=1
				else:
					matrix[row][col] = -1.0 # true negative -> -1
					ntrueneg+=1

	matrix = convert_aa_names(matrix, options.aacode)

	# Generate header with statistics
	npoints=ntrueneg+nfalseneg+ntruepos+nfalsepos
	header = '# Total number of points            : {:6d}\n'.format(npoints+nignored)
	header += '# Number of points ignored   (0.0) : {:6d}\n'.format(nignored)
	header += '# Number of points considered      : {:6d}\n'.format(npoints)
	header += '# Number of true positives   (1.0) : {:6d} ({:3.1f}%)\n'.format(ntruepos, float(ntruepos)/float(npoints)*100.0)
	header += '# Number of false positives  (0.5) : {:6d} ({:3.1f}%)\n'.format(nfalsepos, float(nfalsepos)/float(npoints)*100.0)
	header += '# Number of true negatives  (-1.0) : {:6d} ({:3.1f}%)\n'.format(ntrueneg, float(ntrueneg)/float(npoints)*100.0)
	header += '# Number of false negatives (-0.5) : {:6d} ({:3.1f}%)\n'.format(nfalseneg, float(nfalseneg)/float(npoints)*100.0)

	# Print out comparison map
	global float_digits
	float_digits = 1
	print_map(matrix, header, options.outfile)


# roc_curve_analysis subcommand --------------------------------------------------
# @brief    compares two contact maps and generates ROC parameters
# @detail   compares the contact map of a structure or set of structures to the
#           contact map of a single reference structure and calculets the
#           parameters for ROC analysis based on frequency of the contact in the
#           set of models
def roc_curve_analysis(options):

	matrix, models, comments = read_contact_map_file(options.native_file, True)
	matrix_mod, models, comments = read_contact_map_file(options.model_file, True)

	# Set PDB-ID string
	pdbid="\n"
	if options.pdbid != "":
		pdbid = '\t'+ options.pdbid + "\n"

	# set Header string
	header ='cutoff\t'  # cutoff Value
	header+='co_perc\t' # cutoff Value (percentage)
	header+='TP\t'      # number of true positives
	header+='FP\t'      # number of false posiives
	header+='TN\t'      # number of true negative
	header+='FN\t'      # number of false negatives
	header+='FPR\t'     # false positive rate
	header+='TPR'       # true positive rate
	if options.statistics :
		header+='\tACC\t'   # accuracy
		header+='SPC\t'     # specificity
		header+='PPV\t'     # positive predictive value
		header+='NPV\t'     # negative predictive value
		header+='FDR\t'     # false discovery rate
		header+='MCC\t'     # Matthews correlation coefficient
		header+='F1S'       # F1 score
	if options.pdbid != "":
		header+='\tpdb'
	header +='\n'

	# Get the highest contact count of the model map and scale cutoffs accordingly
	max_val =float( max((max(x[1:], key=float) for x in matrix_mod[1:]), key=float))
	if max_val > 1:
		max_val = (int(max_val) / 10 +1) *10
	else:
		max_val=1
	cutoffs = [x / 100.0 * float(max_val) for x in range(0,101)]

	#Set up parameters for the area under curve calculation
	trueposrate=0.0
	falseposrate=0.0
	auc=0.0
	output_string=""

	# Calculate parameters for every cutoff value
	for cutoff in cutoffs :
		nrow=0
		nignored=0
		ntruepos=0
		nfalsepos=0
		ntrueneg= 0
		nfalseneg= 0

		# Loop over the fields in the contact matrices
		for row in range(1,min(len(matrix), len(matrix_mod))):
			for col in range(1,min(len(matrix[1]), len(matrix_mod[1]))):
				# Skip entries closer in sequence than the specified sequence separation value
				if (math.fabs(col-row)<options.seqsep ):
					nignored +=1
				#
				elif(options.exclude and matrix_mod[row][col]==0):
					nignored +=1
				# Check if value in model matrix is considered positive based on the current cutoff
				elif (float(matrix_mod[row][col])>=cutoff):
					# Compare to the reference structure and increase the respective count
					if (matrix[row][col] == 1):
						ntruepos += 1 # true positive
					else:
						nfalsepos += 1 # false positive
				else:
					if (matrix[row][col] == 1):
						nfalseneg+=1 # false negative
					else:
						ntrueneg+=1 # true negative

		# Calculate statistical parameters for the current cutoff and add them to the output
		output_string+='{}\t{:3.2f}\t{}\t{}\t{}\t{}'.format(cutoff, float(cutoff)/max_val,ntruepos,nfalsepos,ntrueneg,nfalseneg)
		# true positive rate
		tpr=0
		if ntruepos+nfalseneg >0 :
			tpr = float(ntruepos)/(ntruepos+nfalseneg)
		# false positive rate
		fpr = float(nfalsepos)/(nfalsepos+ntrueneg)
		# add current increment to the area under the curve (auc)
		if (cutoff > 0):
			auc = auc + 0.5 * (tpr+trueposrate)*(falseposrate-fpr)
		# Store rates for calculation of next auc step
		trueposrate=tpr
		falseposrate=fpr
		output_string+='\t{:4.3f}\t{:4.3f}'.format(fpr, tpr)
		if options.statistics :
			# accuracy + specificity
			npoints=ntrueneg+nfalseneg+ntruepos+nfalsepos
			accuracy=float(ntruepos+ntrueneg)/npoints
			specificity=1-fpr
			# positive predictive value
			pospredval=0
			if (ntruepos+nfalsepos>0):
				pospredval=ntruepos/float(ntruepos+nfalsepos)
			# negative predictive value
			negpredval=0
			if(ntrueneg+nfalseneg>0): negpredval=ntrueneg/float(ntrueneg+nfalseneg)
			# false discovery rate
			falsediscrate=0
			if(ntruepos+nfalsepos>0):
				falsediscrate=nfalsepos/float(ntruepos+nfalsepos)
			# Matthews correlation coefficient
			matcorcoeff=0
			if((ntruepos+nfalseneg)*(nfalsepos+ntrueneg)*(ntruepos+nfalsepos)*(ntrueneg+nfalseneg)>0):
				matcorcoeff= float(ntruepos*ntrueneg-nfalsepos*nfalseneg)/math.sqrt((ntruepos+nfalseneg)*(nfalsepos+ntrueneg)*(ntruepos+nfalsepos)*(ntrueneg+nfalseneg))
			# F1 score
			fscore=0
			if ntruepos+nfalsepos+nfalseneg>0:
				fscore=2*ntruepos/float(2*ntruepos+nfalsepos+nfalseneg)
			# Add statistics to output
			output_string += '\t{:3.2f}\t{:3.2f}\t{:3.2f}\t{:3.2f}\t{:3.2f}\t{:3.2f}\t{:3.2f}'.format(accuracy,specificity,pospredval,negpredval,falsediscrate,matcorcoeff,fscore)
		output_string += pdbid
	fh = sys.stdout if options.outfile =='' else open(options.outfile, 'w')
	fh.write('# AUC:\t{:3.2f}\n'.format(auc))
	fh.write(header+output_string)
	fh.flush()
	if options.outfile !='':
		fh.close()


# Renumber Subcommand ----------------------------------------------------------
# @brief    Renumbers the headers contact map
# @detail   Renumbers the header row and column of a contact map
#           (default matrix format) from Rosetta numbering to a more appropriate
#           numbering and/or converts the default three-letter amino acid code
#           to one letter code if desired
def renumber(options):

	# Initialize Variables
	row_start, col_start = parse_option_string(options.start)
	row_offset, col_offset = parse_option_string(options.offset)

	# Loop over input files
	for curr_file in options.input_files :
		matrix, models, comments = read_contact_map_file(curr_file)

		# Process header row if specified
		if options.pos == 'col' or options.pos == 'both' :
			col_pos = col_start
			# Loop over entries
			for col in range (1, len(matrix[0])):
				# Split string into residue name and residue number and process them accordingly
				residue_string = matrix[0][col]
				split_pos=re.search('[0-9]',residue_string).start(0)
				residue_name =residue_string[:split_pos]
				residue_id=int(residue_string[split_pos:])
				if residue_id >= col_offset:
					matrix[0][col] = residue_name+str(col_pos)
					col_pos+=1

		# Process header column if specified
		if options.pos == 'row' or options.pos == 'both' :
			row_pos = row_start
			for row in range(1, len (matrix)):
				# Split string into residue name and residue number and process them accordingly
				residue_string = matrix[row][0]
				split_pos=re.search('[0-9]',residue_string).start(0)
				residue_name = residue_string[:split_pos]
				residue_id=int(residue_string[split_pos:])
				if residue_id >= row_offset:
					matrix[row][0] = residue_name+str(row_pos)
					row_pos+=1


		matrix = convert_aa_names(matrix, options.aacode)

		# Print the renumbered file
		header='# Number of Models: {}'.format(models)
		if comments != '':
			header += '\n' + comments
		print_map(matrix, header, curr_file)


# Convert Subcommand -----------------------------------------------------------
# @brief    Convert the amino acid code of the headers of a contact map
# @detail   Converts the default three-letter amino acid code to one letter code
#           and vice versa
def convert(options):

	# Loop over input files
	for curr_file in options.input_files :
		matrix, models, comments = read_contact_map_file(curr_file)
		# Set default aa conversion for convert subcommand to 1
		if options.aacode == 0:
			options.aacode = 1

		# Rename Residues
		matrix = convert_aa_names(matrix, options.aacode)

		# Print the modified file
		header='# Number of Models: {}'.format(models)
		if comments != '':
			header += '\n' + comments
		print_map(matrix, header, curr_file)


# Call main function -----------------------------------------------------------
if __name__ == "__main__":
	main()



