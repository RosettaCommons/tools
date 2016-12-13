import sys
import os
from string import join

# Move to utils.py
def log(msg):
	print(msg)

def prune_leading_lines(lines):
	lines = lines[BOILERPLATE_LENGTH:]
	start = 0
	for l in lines:
		sl = l.strip()
		if sl and sl[0:2] != "//" and sl[0:2] != "/*":
			break
		start += 1
	return lines[start:]

def prune_trailing_lines(lines):
	idx = len(lines) - 1
	while idx >= 0:
		sl = lines[idx].strip()
		if sl and sl[0:2] != "//" and sl[0:2] != "/*":
			break
		idx -=1
	return lines[:idx + 1]

def prune(filedata):
	if not filedata.get("prunedlines"):
		filedata["prunedlines"] = prune_trailing_lines(prune_leading_lines(filedata["lines"]))
	return filedata["prunedlines"]

def getAllFiles(folder):
	listOfFiles = []
	filelist = os.listdir(folder)
	for f in filelist:
		fullpath = os.path.join(folder, f)
		if os.path.isdir(fullpath) and f != ".svn":
			listOfFiles.extend(getAllFiles(fullpath))
		elif os.path.isfile(fullpath):
			listOfFiles.append(fullpath)
	return listOfFiles

def getListOfFiles(root, extension_masks = None):
	'''Recursively get a list of all files under root then filter on extension_masks'''
	pathjoin = os.path.join
	filelist = []

	log("Finding files.")

	# print a # or some other character to indicate progress every nth directory/file
	filelist = getAllFiles(root)

	#todo: filter on extensions

	return filelist

# Move to constants.py
BOILERPLATE = [  '''// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-''',
			'''// vi: set ts=2 noet:''',
			'''//''',
			'''// (c) Copyright Rosetta Commons Member Institutions.''',
			'''// (c) This file is part of the Rosetta software suite and is made available under license.''',
			'''// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.''',
			'''// (c) For more information, see http://www.rosettacommons.org. Questions about this can be''',
			'''// (c) addressed to University of Washington CoMotion, email: license@uw.edu.''',
		]
BOILERPLATE_LENGTH = len(BOILERPLATE)

PYBOILERPLATE = [  '''#!/usr/bin/env python''',
			'''#''',
			'''# (c) Copyright Rosetta Commons Member Institutions.''',
			'''# (c) This file is part of the Rosetta software suite and is made available under license.''',
			'''# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.''',
			'''# (c) For more information, see http://www.rosettacommons.org. Questions about this can be''',
			'''# (c) addressed to University of Washington CoMotion, email: license@uw.edu.''',
		]
PYBOILERPLATE_LENGTH = len(PYBOILERPLATE)

BannedWords = [
	"goto",
	"scanf",
]

def checkHeader(filedata):
	index = 0
	lines = filedata["lines"]
	for hline in BOILERPLATE:
		if lines[index] != BOILERPLATE[index]:
			return["The file header does not match the standard header (emacs settings and copyright information):\n  %s\n  %s" % (lines[index], BOILERPLATE[index])]
		index += 1
	return []

def checkPyHeader(filedata):
	index = 0
	lines = filedata["lines"]
	for hline in PYBOILERPLATE:
		if lines[index] != PYBOILERPLATE[index]:
			return["The file header does not match the standard header (settings and copyright information):\n  %s\n  %s" % (lines[index], PYBOILERPLATE[index])]
		index += 1
	return []

def checkBannedWords(filedata):
	warnings = []
	words = filedata["words"]
	found = [b for b in BannedWords if b in words]
	if found:
		warnings.append("Banned word(s) found: %s." % (join(found, ",")))
	return warnings

def checkGuards(filedata):
	prune(filedata)
	lines = filedata["prunedlines"]
        if len(lines) == 0:
            #No parsed lines means a trivial header that doesn't need includes (e.g. documentation only)
            return []
	idx = 0
	for l in lines:
		if l and not(l.startswith("#include")):
			break
		idx += 1
        if idx >= len(lines):
            #A "meta" include - only has include lines,. The sub includes should have their own include guards.
            return []
	if not(lines[idx].startswith("#ifndef") and lines[-1].startswith("#endif")):
		return ["Missing inclusion guards."]
	return []

def getIncGuardName(filename):
	path, name = os.path.split(filename)
	while path:
		path, curdir = os.path.split(path)
		if curdir == 'src':
			break
		name = curdir + '_' + name
	if not path:
		raise ValueError("Directory 'src' not found in path " + filename)
	name = "INCLUDED_" + name
	name = name.replace('.','_')
	return name

def checkGuardsWithNaming(filedata):
	prune(filedata)
	lines = filedata["prunedlines"] # Should have blank and comment lines removed from front and back
	if len(lines) <= 2:
		return ["Abnormally short file!"]
	if not lines[0].startswith("#ifndef"):
		return ["Proper include guard not found"]
	if not lines[-1].startswith("#endif"):
		return ["Include guard endif not found"]
	ifndef_name = lines[0].split()[1]
	# Capitalization not really important (shouldn't have cap-dependent paths).
	if ifndef_name.lower() != getIncGuardName(filedata["filename"]).lower():
		return ["Include guard not named appropriately ("+ifndef_name+")."]
	if not lines[1].startswith("#define"):
		return ["No define for include guard."]
	if lines[1].split()[1] != ifndef_name: # Capitalization *here* is important.
		return ["Include guard define not named appropriately ("+lines[1].split()[1]+")."]
	return []

CodingConventions = [
        (checkHeader, ["hh"]),
        (checkHeader, ["cc"]),
#        (checkPyHeader, ["py"]),
#        (checkBannedWords, ["hh"]),
#        (checkGuards, ["hh"]),
#        (checkGuardWithNaming, ["hh"]),
    ]

def check_source(root):
	warnings = []
	filelist = getListOfFiles(root)
	for fl in filelist:
		F = open(fl, "r")
		contents = F.read()
		lines = contents.split("\n")
		words = set(contents.split())
		F.close()
		cwarnings = []

		filedata = {"contents" : contents, "lines" : lines, "words" : words, "filename" : os.path.abspath(fl) }
		log("Processing %s:" % fl)
		for convention in CodingConventions:
			for ext in convention[1]: #Think of a better solution here
				if fl.endswith(ext):
					cwarnings.extend(convention[0](filedata))
					break
		if cwarnings:
			warnings.append((fl, cwarnings))

	return warnings

def burninate(root):
	print(root)
	warnings = check_source(root)
	if warnings:
		for file, wrngs in warnings:
			log("%s:" % file)
			for w in wrngs:
				log("\t%s" % w)

	log("Done")

if __name__ == "__main__":
	if len(sys.argv) == 1:
		burninate(os.path.join(os.getcwd(), "src", "numeric"))
	elif len(sys.argv) == 2 and sys.argv[1] != "-h":
		burninate(os.path.abspath(sys.argv[1]))
	else:
		print "Usage:", sys.argv[0], " [ path ] "
		print "(Path defaults to ./src/numeric)"
