#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import urllib2 as url
# The Bio module is provided by BioPython
from Bio.PDB import PDBParser as pdbp
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB import PDBIO
import re

if len(sys.argv) < 3:
	print "convert_pdb_to_antibody_numbering_scheme.py pdb output_pdb <heavy_chain> <light_chain> <scheme>\n \
	ex. convert_pdb_to_antibody_numbering_scheme.py 3NGB.pdb 3NGB_kabat.py H L c\n \
	sheme = c for chothia, k for kabat, e for extended chothia\n"
	sys.exit()



pdbparser = pdbp()
structure = pdbparser.get_structure('self',sys.argv[1])
output = sys.argv[2]


try:
	h_chain = sys.argv[3]

except IndexError:
	h_chain = "H"


try:
	l_chain = sys.argv[4]
except IndexError:
	l_chain = "L"

heavy_chain_res = structure[0][h_chain]
light_chain_res = structure[0][l_chain]

try:
	scheme = sys.argv[5]
except IndexError:
	scheme = "c"


def get_sequence(chain_object):
	res_names = []
	for residue in chain_object:
		res = three_to_one(residue.resname)
		res_names.append(res)
	return "".join(res_names)

def get_numbering_scheme(ab_string,scheme=scheme):
	url_ = "http://www.bioinf.org.uk/abs/abnum/abnum.cgi?plain=1&aaseq={}&scheme=-{}".format(ab_string,scheme)
	response = url.urlopen(url_)
	scheme = []
	for i in response.readlines():
		scheme_number = i.strip().split()[0][1:]
		matcher = i.strip().split()[1]
		if matcher == "-":
			continue
		insertion_code = ' '
		re_ = re.search('[A-Za-z]',scheme_number)
		if re_:
			insertion_code = re_.group(0)
			scheme_number = scheme_number[:scheme_number.index(insertion_code)]
		scheme_number = int(scheme_number)
		scheme.append((' ',scheme_number,insertion_code))
	return scheme

def set_new_res(resis):
	tmp_dict = {}
	for res_,id_ in zip(resis,get_numbering_scheme(get_sequence(resis),scheme=scheme)):
		tmp_dict[res_.id] = id_
		resnum = int(id_[1]) #find where the url stops coding at, then we pick it up from there
	for res_ in resis:
		try:
			res_.id = tmp_dict[res_.id]
		except KeyError:
			resnum += 1
			res_.id = ( ' ',resnum, ' ')


set_new_res(heavy_chain_res)
set_new_res(light_chain_res)

w = PDBIO()
w.set_structure(structure)
w.save(output)
