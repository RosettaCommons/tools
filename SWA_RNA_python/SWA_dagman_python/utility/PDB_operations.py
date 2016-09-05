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
#######################


from SWA_util import *


####################################################################################
###Nov 11, 2011
standard_polar_hatom_list=["HO2'", " H21", " H22", " H41", " H42", " H61", " H62", " H1 ", " H3 ", "HO5'", "HO3'"] #"HO5'", "HO3'" are (Non-Rosetta) variants of 3' and 5' chain-ends

standard_nonpolar_hatom_list=[" H1'", " H2'", " H3'", " H4'", " H5'", "H5''", " H2 ", " H5 ", " H6 ", " H8 "]

###########################################################################################
###Nov 11, 2011
def Is_standard_polar_hatom_name(query_atom_name, matching_whitespace):

	if(matching_whitespace):
		
		if(query_atom_name in standard_polar_hatom_list):
			return True
		else:
			return False

	else:

		for hatom_name in standard_polar_hatom_list:

			if(len(query_atom_name.split())!=1): error_exit_with_message("len(query_atom_name.split())!=1")
			if(len(hatom_name.split())!=1): error_exit_with_message("len(hatom_name.split())!=1")

			if(query_atom_name.split()[0]==hatom_name.split()[0]): return True

		return False

###########################################################################################
###Nov 11, 2011
def Is_standard_nonpolar_hatom_name(query_atom_name, matching_whitespace):

	if(matching_whitespace):
		
		if(query_atom_name in standard_nonpolar_hatom_list):
			return True
		else:
			return False

	else:

		for hatom_name in standard_nonpolar_hatom_list:

			if(len(query_atom_name.split())!=1): error_exit_with_message("len(query_atom_name.split())!=1")
			if(len(hatom_name.split())!=1): error_exit_with_message("len(hatom_name.split())!=1")

			if(query_atom_name.split()[0]==hatom_name.split()[0]): return True

		return False

###########################################################################################
###Nov 11, 2011
def Is_standard_hatom_name(query_atom_name, matching_whitespace):

	if(Is_standard_polar_hatom_name(query_atom_name, matching_whitespace)): return True

	if(Is_standard_nonpolar_hatom_name(query_atom_name, matching_whitespace)): return True
	
	return False

###########################################################################################

###Nov 11, 2011
def assert_is_valid_standard_hatom_name(query_atom_name, matching_whitespace):

	if(Is_standard_hatom_name(query_atom_name, matching_whitespace)==False): error_exit_with_message("query_atom_name (%s) is not a valid standard_hatom_name!" %(query_atom_name))
###########################################################################################



######Oct 14, 2011, WARNING THIS IS INDEPENDENT of replace_pdb_in_place.py##########
def renumber_pdb_in_place_func(input_pdb):

	if(exists(input_pdb)==False): error_exit_with_message("input_pdb (%s) doesn't exist!")

	lines = open(input_pdb,'r').readlines()

	out = open(input_pdb,'w') #Overwrite the old pdb_file!

	oldresnum = '   '
	oldchainID = ' '
	oldlongname = '  '
	count = 0;

	outid  = open( input_pdb,'w') #Overwrite the old pdb_file!

	atomnum  = 0
	for line in lines:
		line_edit = line
		if(line[0:3] == 'TER'): continue

		if(line_edit.count('HETATM')>0): error_exit_with_message("line_edit.count('HETATM')>0")

		if(line.count('ATOM')>0):

			if(line[0:4] != 'ATOM'): error_exit_with_message("line.count('ATOM')>0 but line[0:4] != 'ATOM'")

			##line[16] is the Alternative location indicator, check that there is only one conformation (i.e. 'A' or empty!)
			if(line[16]!=' ' and line[16]!='A'): error_exit_with_message("line[16]!=' ' and line[16]!='A', line[16]=%s" %(line[16]))

			atomnum += 1

			resnum = line_edit[22:26] #Change from 23:26 to 22:26 on Nov 11, 2011
			chainID=line_edit[21]
			longname = line_edit[17:20]
			insertcode=line_edit[26]

			if(insertcode!= ' '): error_exit_with_message("insertcode!= ' ' for line=%s" %(line))

			if( (resnum != oldresnum) or (chainID!=oldchainID) ): 
				count = count + 1
			else:
				if(longname!=oldlongname): error_exit_with_message("longname!=oldlongname")
				

			oldresnum  = resnum
			oldchainID = chainID
			oldlongname= longname

			newnum = '%4d' % count
			line_edit = '%s%5d%s%s%s' % (line_edit[0:6], atomnum, line[11:22], newnum, line_edit[26:] )

			outid.write(line_edit) 

	outid.close()


######Oct 14, 2011, WARNING THIS IS INDEPENDENT of replace_chain_inplace.py########
def replace_chain_inplace_func(input_pdb, newchain):

	if(len(newchain)!=1): error_exit_with_message("len(newchain)!=1, newchain=%s" %(newchain))

	if(exists(input_pdb)==False): error_exit_with_message("input_pdb (%s) doesn't exist!")

	lines = open(input_pdb,'r').readlines()

	out = open(input_pdb,'w') #Overwrite the old pdb_file!

	for i in range( len(lines)):
		line = lines[i]
		if(line.count('ATOM')>0):

			if(line[0:4] != 'ATOM'): error_exit_with_message("line.count('ATOM')>0 but line[0:4] != 'ATOM'")

			line = line[0:21]+newchain+line[22:]
			out.write(line)
		else: #New Oct 14, 2011..just print out the old_line
			out.write(line)



############################################################################################################


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

		################################Oct 19, 2011###################################################################
		elif( (line[19]=="A") and (line[12:16]=="H1  " or line[12:16]==" H1 ") ): #Special case, protonated A!
			#"H1  " is for backward compatibility due to a error!
			new_line_list[77]="H"
			new_line_list[12:16]=" H1 "
		################################################################################################################

		################################Variants of 3' and 5' chain-ends################################################
		elif(line[12:16]=="HO5*"): #,"5HO*"
			new_line_list[77]="H"
			new_line_list[12:16]="HO5'"
		elif(line[12:16]=="HO3*"): #,"3HO*"
			new_line_list[77]="H"
			new_line_list[12:16]="HO3'"
			
		else:
			error_exit_with_message("invalid atom_name=%s, line=%s" %(line[12:16], line) )
	

		#Stop using this on Oct 02, 2011.
		#new_line_list[17]=line[19]
		#new_line_list[18]=' '
		#new_line_list[19]=' '


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



