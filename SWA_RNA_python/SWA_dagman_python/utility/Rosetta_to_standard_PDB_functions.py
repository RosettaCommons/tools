#!/usr/bin/env python

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
import os
from os.path import basename, dirname, exists, expanduser
from time import sleep
import copy
import string
from os import popen 
###########################

from SWA_util import *
###########################


def Rosetta_to_standard_PDB_func(input_pdb_file, remove_hydrogen, output_pdb_file, VERBOSE=True):

	if(input_pdb_file==""): error_exit_with_message("input_pdb_file==\"\"")  

	if(exists(input_pdb_file)==False): error_exit_with_message("input_pdb_file (%s) doesn't exist!!" %(input_pdb_file))  

	#########################################################################################################
	if(output_pdb_file==""): 
		output_pdb_file="STANDARD_%s" %(basename(input_pdb_file))
	else:
		"User inputted output_pdb_file= %s " %(output_pdb_file)
		if(exists(output_pdb_file)): error_exit_with_message("User inputted output_pdb_file (%s) already exist!" %(output_pdb_file) )

	if(exists(output_pdb_file)):
		print "Warning output_pdb_file (%s) already exist...removing!" %(output_pdb_file)
		submit_subprocess( "rm %s " %(output_pdb_file)  )
	#########################################################################################################

	#if(remove_hydrogen==False): output_pdb_file="WITH_H_" + output_pdb_file ###Commented Out on Oct 2, 2011, after chance DEFAULT value of remove_hydrogen to false

	if(remove_hydrogen): output_pdb_file="No_H_" + output_pdb_file ###Oct 2, 2011, after chance DEFAULT value of remove_hydrogen to false


	no_hydrogen_atom_pdb=append_to_basename("No_H_" , input_pdb_file) 

	if(remove_hydrogen):
		submit_subprocess("remove_hydrogen_from_pdb.py -s %s -output_pdb %s" %(input_pdb_file, no_hydrogen_atom_pdb) ) 

		if(exists(no_hydrogen_atom_pdb)==False): error_exit_with_message("no_hydrogen_atom_pdb (%s) doesn't exist!")

		job_string_list=open(no_hydrogen_atom_pdb).readlines()

	else:

		job_string_list=open(input_pdb_file).readlines()



	OUPUT_PDB = open( output_pdb_file, 'w')

	for line in job_string_list:

		if(line=="END\n"): continue 

		if(line[0:3]=="TER"): continue 

		if(line[0:6]!="ATOM  "): error_exit_with_message("Not a ATOM line, line=%s" %(line))
		
		#if(len(line)!= 81): error_exit_with_message("len(line)!= 81, len(line)=%d, line=%s" %(len(line), line) )  


		if(line[18]!='r'): error_exit_with_message("line[18]!='r', line=%s" %(line) )  
		if(line[19]!='G' and line[19]!='U' and line[19]!='C' and line[19]!='A'): error_exit_with_message("line[19]!='G' and line[19]!='U' and line[19]!='C' and line[19]!='A' line=%s" %(line))
		if(line[20]!=' '): error_exit_with_message("line[20]!=' ', line=%s" %(line) )  

		if(remove_hydrogen):
			if(line[12:16].count("H")>0): error_exit_with_message("remove_hydrogen mode but encounter hydrogen atom (%s) in line=%s" %(line[12:16], line) )

		new_line_list=[]

		for char in line:	
			new_line_list.append(char)

		if(line[12:16]==" P  "): new_line_list[77]="P"

		elif(line[12:16]==" C1*"): new_line_list[77]="C"
		elif(line[12:16]==" C2*"): new_line_list[77]="C"
		elif(line[12:16]==" C3*"): new_line_list[77]="C"
		elif(line[12:16]==" C4*"): new_line_list[77]="C"
		elif(line[12:16]==" C5*"): new_line_list[77]="C"

		elif(line[12:16]==" C2 "): new_line_list[77]="C"
		elif(line[12:16]==" C4 "): new_line_list[77]="C"
		elif(line[12:16]==" C5 "): new_line_list[77]="C"
		elif(line[12:16]==" C6 "): new_line_list[77]="C"
		elif(line[12:16]==" C8 "): new_line_list[77]="C"

		elif(line[12:16]==" O1P"): new_line_list[77]="O"
		elif(line[12:16]==" O2P"): new_line_list[77]="O"
		elif(line[12:16]==" O2*"): new_line_list[77]="O"
		elif(line[12:16]==" O3*"): new_line_list[77]="O"
		elif(line[12:16]==" O4*"): new_line_list[77]="O"
		elif(line[12:16]==" O5*"): new_line_list[77]="O"

		elif(line[12:16]==" O2 "): new_line_list[77]="O"
		elif(line[12:16]==" O4 "): new_line_list[77]="O"
		elif(line[12:16]==" O6 "): new_line_list[77]="O"

		elif(line[12:16]==" N1 "): new_line_list[77]="N"
		elif(line[12:16]==" N2 "): new_line_list[77]="N"
		elif(line[12:16]==" N3 "): new_line_list[77]="N"
		elif(line[12:16]==" N4 "): new_line_list[77]="N"
		elif(line[12:16]==" N6 "): new_line_list[77]="N"
		elif(line[12:16]==" N7 "): new_line_list[77]="N"
		elif(line[12:16]==" N9 "): new_line_list[77]="N"

		############Oct 02, 2011..consider the hydrogen atoms#########
		elif(line[12:16]=="1H5*"): new_line_list[77]="H"
		elif(line[12:16]=="2H5*"): new_line_list[77]="H"
		elif(line[12:16]==" H4*"): new_line_list[77]="H"
		elif(line[12:16]==" H3*"): new_line_list[77]="H"
		elif(line[12:16]==" H1*"): new_line_list[77]="H"
		elif(line[12:16]=="1H2*"): new_line_list[77]="H"
		elif(line[12:16]=="2HO*"): new_line_list[77]="H"
		elif(line[12:16]==" H1 "): new_line_list[77]="H"  #imino group of G or protonated A
		elif(line[12:16]==" H3 "): new_line_list[77]="H"  #imino group of U or protonated C
		elif(line[12:16]=="1H2 "): new_line_list[77]="H"  #amino group of G
		elif(line[12:16]=="2H2 "): new_line_list[77]="H"  #amino group of G
		elif(line[12:16]=="1H4 "): new_line_list[77]="H"  #amino group of C
		elif(line[12:16]=="2H4 "): new_line_list[77]="H"  #amino group of C
		elif(line[12:16]=="1H6 "): new_line_list[77]="H"  #amino group of A
		elif(line[12:16]=="2H6 "): new_line_list[77]="H"  #amino group of A
		elif(line[12:16]==" H2 "): new_line_list[77]="H"  #C-H group of A
		elif(line[12:16]==" H5 "): new_line_list[77]="H"  #C-H group of C or U
		elif(line[12:16]==" H6 "): new_line_list[77]="H"  #C-H group of C or U
		elif(line[12:16]==" H8 "): new_line_list[77]="H"  #C-H group of G or A

		########OK one other possibility is 3HO* or HO3' hydrogen at the 3' end of each chain (add by Molprobity), BUT this is not part of the standard Rosetta nucleotide####

		else:
			error_exit_with_message("invalid atom_name=%s, line=%s" %(line[12:16], line) )
	
		'''//Stop using this on Oct 02, 2011.
		new_line_list[17]=line[19]
		new_line_list[18]=' '
		new_line_list[19]=' '
		'''

		##Oct 02 2011: Based on recent PDB structure, the nucleotide name is located at position 19. This might lead to some backward compabilility issue!
		new_line_list[17]=' '
		new_line_list[18]=' '
		new_line_list[19]=line[19]


		if(line[12:16]==" O1P" ): new_line_list[12:16]=" OP1"
		elif(line[12:16]==" O2P" ): new_line_list[12:16]=" OP2"
		elif(line[12:16]=="1H5*" ): new_line_list[12:16]=" H5'"
		elif(line[12:16]=="2H5*" ): new_line_list[12:16]="H5''"
		elif(line[12:16]=="1H2*" ): new_line_list[12:16]=" H2'"
		elif(line[12:16]=="2HO*" ): new_line_list[12:16]="HO2'"
		elif(line[12:16]=="1H2 " ): new_line_list[12:16]=" H21"
		elif(line[12:16]=="2H2 " ): new_line_list[12:16]=" H22"
		elif(line[12:16]=="1H4 " ): new_line_list[12:16]=" H41"
		elif(line[12:16]=="2H4 " ): new_line_list[12:16]=" H42"
		elif(line[12:16]=="1H6 " ): new_line_list[12:16]=" H61"
		elif(line[12:16]=="2H6 " ): new_line_list[12:16]=" H62"



		if(line[15]=="*"):  new_line_list[15]="'" #All other atom name with prime w
	
		seperator=""

		new_line=seperator.join(new_line_list)
		#print "new_line= %s" %(new_line)

		if(new_line[12:16].count("*")>0): error_exit_with_message("atom_name contain the * character, for new_line=%s: "%(new_line[12:16], new_line) )

		OUPUT_PDB.write(new_line)

	OUPUT_PDB.close()

	if(remove_hydrogen): submit_subprocess("rm %s " %(no_hydrogen_atom_pdb) )





########ATOM     41  O4*  rG A   2     -10.469  -6.154   2.097  1.00 10.00              
########ATOM     24 1H5*  rG A   1     -12.562  -6.584   9.293  1.00  0.00  
########12345678901234567890123456789012345678901234567890123456789012345678901234567890
########0        1         2         3         4         5         6         7         8





'''
01234567890123456789012345678901234567890123456789012345678901234567890123456789	
ATOM      1  P    rG R   1     101.133 -43.649  16.700  1.00 59.94           P  
ATOM      2  O1P  rG R   1     101.879 -44.726  16.006  1.00 60.94           O  
ATOM      3  O2P  rG R   1     101.834 -42.388  17.023  1.00 56.97           O  
ATOM      4  O5*  rG R   1      99.855 -43.236  15.828  1.00 56.83           O  
ATOM      5  C5*  rG R   1      99.271 -44.139  14.899  1.00 55.00           C  
ATOM      6  C4*  rG R   1      98.207 -43.454  14.054  1.00 53.51           C  
ATOM      7  O4*  rG R   1      97.351 -42.626  14.888  1.00 51.34           O  
ATOM      8  C3*  rG R   1      98.731 -42.485  13.000  1.00 52.08           C  
ATOM      9  O3*  rG R   1      99.328 -43.159  11.846  1.00 51.13           O  
ATOM     10  C2*  rG R   1      97.455 -41.706  12.689  1.00 51.48           C  
ATOM     11  O2*  rG R   1      96.588 -42.362  11.784  1.00 51.16           O  
ATOM     12  C1*  rG R   1      96.835 -41.578  14.089  1.00 52.18           C  
ATOM     13  N9   rG R   1      97.175 -40.288  14.674  1.00 52.27           N  
ATOM     14  C8   rG R   1      98.138 -40.006  15.612  1.00 52.60           C  
ATOM     15  N7   rG R   1      98.207 -38.733  15.905  1.00 54.08           N  
ATOM     16  C5   rG R   1      97.239 -38.137  15.109  1.00 51.50           C  
ATOM     17  C6   rG R   1      96.861 -36.780  14.975  1.00 51.23           C  
ATOM     18  O6   rG R   1      97.309 -35.795  15.575  1.00 49.84           O  
ATOM     19  N1   rG R   1      95.826 -36.601  14.053  1.00 51.01           N  
ATOM     20  C2   rG R   1      95.261 -37.624  13.321  1.00 52.69           C  
ATOM     21  N2   rG R   1      94.284 -37.269  12.470  1.00 51.75           N  
ATOM     22  N3   rG R   1      95.605 -38.902  13.435  1.00 50.38           N  
01234567890123456789012345678901234567890123456789012345678901234567890123456789	
0         1         2         3         4         5         6         7  

#Convert 1 to 2:
#					ATOM      2  C5*  rG A   1      40.311  50.542  48.455  1.00 33.79           C  (Rosetta)
#					ATOM      2  C5' G   B   1      40.311  50.542  48.455  1.00 33.79           C  (ModeRNA)

#					ATOM     22  O1P  rG A   1      41.587  49.843  50.847  1.00 43.84           O  (Rosetta)
#					ATOM     22  OP1 G   B   1      41.587  49.843  50.847  1.00 43.84           O  (ModeRNA)

#		   		01234567890123456789012345678901234567890123456789012345678901234567890123456789	
#					0         1         2         3         4         5         6         7     
'''
