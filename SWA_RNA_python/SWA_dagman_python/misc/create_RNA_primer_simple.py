#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
import os
import copy
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################



primer_length=60

RNA_map={}

T7_promoter='TTCTAATACGACTCACTATA'

RNA_map=[]

RNA_map.append(["Herm_control", T7_promoter +'GGAACAAACTTGACAggagTggccgaaaggcaTcTccAGCTAAAAGAAACAACAACAACAAC'])

RNA_map.append(["AA_5'of_loop", T7_promoter +'GGAACAAACUUGACAuacgaagccgaaaggcucguaAGCUAAAAGAAACAACAACAACAAC']) 
RNA_map.append(["UU_5'of_loop", T7_promoter +'GGAACAAACUUGACAuacguugccuucgggcacguaAGCUAAAAGAAACAACAACAACAAC']) 
RNA_map.append(["AA_3'of_loop", T7_promoter +'GGAACAAACUUGACAuacgugccgaaaggcaacguaAGCUAAAAGAAACAACAACAACAAC']) 
RNA_map.append(["UU_3'of_loop", T7_promoter +'GGAACAAACUUGACAuacgagccuucgggcuucguaAGCUAAAAGAAACAACAACAACAAC']) 


RNA_map.append(["U1_5'of_loop", T7_promoter +'GGAACAAACUUGACAuacguggccuucgggcccguaAGCUAAAAGAAACAACAACAACAAC'])  
RNA_map.append(["U2_5'of_loop", T7_promoter +'GGAACAAACUUGACAuacggugccuucgggcccguaAGCUAAAAGAAACAACAACAACAAC']) 
RNA_map.append(["No_bulge_one", T7_promoter +'GGAACAAACUUGACAuacguugccuucgggcaacguaAGCUAAAAGAAACAACAACAACAAC'])
RNA_map.append(["No_bulge_two", T7_promoter +'GGAACAAACUUGACAuacgaagccuucgggcuucguaAGCUAAAAGAAACAACAACAACAAC'])



#RNA_map["AA_5'of_loop"]=T7_promoter +'GGAACAAACUUGACAggcgaagccgaaaggcucgccAGCUAAAAGAAACAACAACAACAAC'
#RNA_map["UU_5'of_loop"]=T7_promoter +'GGAACAAACUUGACAggcguugccuucgggcacgccAGCUAAAAGAAACAACAACAACAAC'
#RNA_map["AA_3'of_loop"]=T7_promoter +'GGAACAAACUUGACAggcgugccgaaaggcaacgccAGCUAAAAGAAACAACAACAACAAC'
#RNA_map["UU_3'of_loop"]=T7_promoter +'GGAACAAACUUGACAggcgagccuucgggcuucgccAGCUAAAAGAAACAACAACAACAAC'


for RNA_pair in RNA_map:

	RNA_name=RNA_pair[0]

	RNA_seq=RNA_pair[1]

	forward_primer=RNA_seq[0:60].upper().replace("U", "T") # This is 5' to 3'

	backward_prime_complement=RNA_seq[::-1][0:60].upper().replace("U", "T") #This is 3' to 5'
	backward_primer=""

	for com_base in backward_prime_complement:
		base=""

		if(com_base=="A"): 
			base="T"
		elif(com_base=="T"):
			base="A"
		elif(com_base=="G"):
			base="C"
		elif(com_base=="C"):
			base="G"
		else:
			error_exit_with_message("Invalid com_base=%s "%(com_base) )

		backward_primer+=base	


	print "--------------------------------------------------------------------------------------------"
	print "RNA_name= " , RNA_name
	print "RNA_seq= " , RNA_seq
	print "forward_primer= " , forward_primer
	print "backwar_primer= " , backward_primer

	print "--------------------------------------------------------------------------------------------"


