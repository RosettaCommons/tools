#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
import os
import copy
######################################################################

from SWA_util import *
######################################################################


def get_one_letter_code(three_letter_code):

	if(three_letter_code=="ADE"): return "A"

	if(three_letter_code=="CYT"): return "C"

	if(three_letter_code=="URA"): return "U"

	if(three_letter_code=="URI"): return "U" #Oct 10, 2011 For parsing NUCHEMICS outfile.

	if(three_letter_code=="GUA"): return "G"

	error_exit_with_message("Invalid three_letter_code (%s) "(three_letter_code) )


def get_WC_complement(base):

	if(base=="u"): return "a"
	
	if(base=="a"): return "u"
	
	if(base=="c"): return "g"
	
	if(base=="g"): return "c"
	
	error_exit_with_message("invalid base %s " %(base) )

def Is_perfect_WC_helix(strand_1_sequence_passin,strand_2_sequence_passin):

	strand_1_sequence=strand_1_sequence_passin.lower()

	strand_2_sequence=strand_2_sequence_passin.lower()

	if(len(strand_1_sequence)!=len(strand_2_sequence)): error_exit_with_message("len(strand_1_sequence)!=len(strand_2_sequence)")

	for n in range(len(strand_1_sequence)):

		if(strand_1_sequence[n]!=get_WC_complement(strand_2_sequence[-(n+1)])):
			print "strand_1_sequence=5'-%s-3' and strand_2_sequence=5'-%s-3' do not form a perfect WC helix!" %(strand_1_sequence, strand_2_sequence) 
			return False


	return True
		


