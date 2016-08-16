#!/usr/bin/env python

import string
from sys import argv,stdout, stderr
from os import popen,system
from os.path import exists,dirname,basename,abspath
######################################################################

from SWA_dagman_python.database.SWA_amino_acids import longer_names
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################


###Oct 01, 2011: Compare to the regular make_rna_rosetta_ready.py, this version correctly rename all the atoms [especially hydrogens] into Rosetta convention.######################
###Hence when the PDB is passed into Rosetta, these hydrogen atom are correctly recognize....previously they weren't recoginize and hence rebuilt using idealized geometry##########
###Should stick will using this version from now on!################################################################################################################################

all_models= parse_options( argv, "all_models", "False" ) #There are usually more than 1 struct for NMR PDB!

output_pdb= parse_options( argv, "output_pdb", "" )

if(exists(output_pdb)):
	print "output_pdb %s exist ...removing" %(output_pdb) 
	submit_subprocess("rm %s " %(output_pdb) )

alter_conform= parse_options( argv, "alter_conform", "A" )

if(alter_conform not in ['A', 'B']): error_exit_with_message("alter_conform not in ['A', 'B'] | alter_conform=%s" %(alter_conform))

'''
ignore_conform = ''

if(alter_conform == 'A'):
	ignore_conform = 'B'
else:
	ignore_conform = 'A'
'''

fastaid = stdout #Nov 22, previously: fastaid=stderr

if(len(argv)!=2): error_exit_with_message("len(argv)!=2, leftover_argv=%s " %(list_to_string(argv) ) )

input_pdb = argv[1]

if(exists(input_pdb)==False): error_exit_with_message("input_pdb (%s) doesn't exist!" %(input_pdb))

#######################Actual script################


print 'Reading ... '+ input_pdb

if(input_pdb[-3:] == '.gz' ): error_exit_with_message("input_pdb[-3:] == '.gz' | input_pdb=%s" %( input_pdb ) )

if(input_pdb[-4:] != '.pdb'): error_exit_with_message("input_pdb[-4:] == '.pdb' | input_pdb=%s" %( input_pdb ) )

lines = open(input_pdb,'r').readlines()

goodnames = [' rA',' rC',' rG',' rU']

new_struct=True
model_count=0

for line in lines:

	if( (len(line)>5) and (line[:6]=='ENDMDL') ):

		if(all_models==True):
			new_struct=True
			continue
		else:
			break

	######################################################################
	if(new_struct==True):
		model_count+=1

		fastaid.write('\n')

		if(output_pdb==""):
			outfile = input_pdb
			outfile = dirname( abspath(outfile) ) + '/' + string.lower( basename(outfile) )

			if(outfile[-4:] != '.pdb'): error_exit_with_message("outfile[-3:] == '.pdb' | outfile=%s" %( outfile ) )

			outfile = outfile[:-4] #Remove .pdb
			outfile += '_RNA_' + alter_conform

			if(all_models): outfile +="_%d" %(model_count)

			outfile +=  '.pdb'
		else:
			if(all_models): error_exit_with_message("all_models==True but user specified output_pdb NAME!")

			outfile = output_pdb

		outid = open( outfile, 'w')

		print 'Writing ... '+ outfile

		fastaid.write('>'+outfile+'\n');

		oldresnum = '   '
		oldchainID = ' '
		oldlongname = '  '
		oldinsertcode= ' '
		count = 0;

		new_struct=False

	#######################################################################

	if(len(line) <= 21):  continue

	#if(line[21] in chainids or ignore_chain): #Remove this on Nov 22, 2011.
	line_edit = line

	if(line[0:3] == 'TER'): continue

	if(line[0:6] == 'HETATM'):

		'''//Comment out on Nov 03, 2011! Since Protien Stuff! 
		if(line[17:20]=='MSE'): #Selenomethionine
			line_edit = 'ATOM  '+line[6:17]+'MET'+line[20:]
			if(line_edit[12:14] == 'SE'):
				line_edit = line_edit[0:12]+' S'+line_edit[14:]
			if((len(line_edit)>75) and line_edit[76:78] == 'SE'):
				line_edit = line_edit[0:76]+' S'+line_edit[78:]
		'''

		if(line[17:20]=='5BU'): #Selenomethionine
			line_edit = 'ATOM  '+line[6:17]+'  U'+line[20:]
		elif(line[17:20]=='OMC'): #Selenomethionine
			line_edit = 'ATOM  '+line[6:17]+'  C'+line[20:]
		elif(line[17:20]=='5MC'):
			line_edit = 'ATOM  '+line[6:17]+'  C'+line[20:]
		elif(line[17:20]==' DC'): #Selenomethionine
			line_edit = 'ATOM  '+line[6:17]+'  C'+line[20:]
		elif(line[17:20]=='CBR'): #Selenomethionine
			line_edit = 'ATOM  '+line[6:17]+'  C'+line[20:]
		elif(line[17:20]=='CB2'):
			line_edit = 'ATOM  '+line[6:17]+'  C'+line[20:]
		elif(line[17:20]=='2MG'):
			line_edit = 'ATOM  '+line[6:17]+'  G'+line[20:]
		elif(line[17:20]=='H2U'):
			line_edit = 'ATOM  '+line[6:17]+'  U'+line[20:]
		elif(line[17:20]=='PSU'):
			line_edit = 'ATOM  '+line[6:17]+'  U'+line[20:]
		elif(line[17:20]=='5MU'):
			line_edit = 'ATOM  '+line[6:17]+'  U'+line[20:]
		elif(line[17:20]=='OMG'):
			line_edit = 'ATOM  '+line[6:17]+'  G'+line[20:]
		elif(line[17:20]=='7MG'):
			line_edit = 'ATOM  '+line[6:17]+'  G'+line[20:]
		elif(line[17:20]=='1MG'):
			line_edit = 'ATOM  '+line[6:17]+'  G'+line[20:]
		elif(line[17:20]==' YG'):
			line_edit = 'ATOM  '+line[6:17]+'  G'+line[20:]
		elif(line[17:20]=='1MA'):
			line_edit = 'ATOM  '+line[6:17]+'  A'+line[20:]
		#####COPY FROM make_rna_rosetta_ready.py on Nov 03, 2011||Aug 15, 2011 Parin additional modified NTS in 23rRNA START
		elif(line[17:20]=='OMU'):
			line_edit = 'ATOM  '+line[6:17]+'  U'+line[20:]
		elif(line[17:20]=='UR3'):
			line_edit = 'ATOM  '+line[6:17]+'  U'+line[20:]
		else:
			print "WARNING... ignoring HETATM line=%s " %(line[:-1])
		#####COPY FROM make_rna_rosetta_ready.py on Nov 03, 2011||Aug 15, 2011 Parin additional modified NTS in 23rRNA END


	if(line_edit[0:4] == 'ATOM'):

		#Don't save alternative conformations.
		if( ( line[16]!=' ') and (line[16] != alter_conform) ): ###Before Nov 21, 2011, used to be (line[16]==ignore_conform)
			print "Alternative conformation line=%s" %(line[:-1])
			continue

		resnum = line_edit[22:26] #Change from 23:26 to 22:26 on Nov 11, 2011
		chainID=line_edit[21]
		insertcode=line_edit[26]

		if( (resnum != oldresnum) or (chainID!= oldchainID) or (insertcode!=oldinsertcode) ):

			count = count + 1
			longname = line_edit[17:20]
			if  (longname in ['  G', 'G  ', 'GUA']): longname = ' rG'
			elif(longname in ['  A', 'A  ', 'ADE']): longname = ' rA'
			elif(longname in ['  C', 'C  ', 'CYT']): longname = ' rC'
			elif(longname in ['  U', 'U  ', 'URA']): longname = ' rU'
			else:
				if(longname not in goodnames): continue

			if longer_names.has_key(longname):
				fastaid.write( longer_names[longname] );
			else:
				error_exit_with_message("invalid longname (%s)" %(longname))
				#fastaid.write( 'X')

		else:
			if(longname!=oldlongname): error_exit_with_message("longname!=oldlongname")

		oldresnum = resnum
		oldchainID= chainID
		oldlongname=longname
		oldinsertcode=insertcode

		if(longname not in goodnames): error_exit_with_message("longname (%s) not in goodnames (%s)" %(longname, list_to_string(goodnames) ) )

		newnum = '%4d' %(count)

		#line_edit = line_edit[0:16] + ' ' + longname + line_edit[20:22] + newnum + line_edit[26:]

		###Nov 22, 2011: Column 27 which is 26 in python numbering is code for insertion of residue########
		###This is now considered as a new nts and hence line_edit[26] should be ' '#######################
		line_edit = line_edit[0:16] + ' ' + longname + line_edit[20:22] + newnum + ' ' + line_edit[27:] 

		line_edit = line_edit.replace('\'','*')

	 	######################################Oct 28, 2011#####################################################

		new_atom_name=""
		if  (line_edit[12:16]==" H5*"): new_atom_name="1H5*"
		elif(line_edit[12:16]=="H5**"): new_atom_name="2H5*"
		elif(line_edit[12:16]==" H2*"): new_atom_name="1H2*"
		elif(line_edit[12:16]=="HO2*"): new_atom_name="2HO*"

		elif(line_edit[12:16]==" H21"): new_atom_name="1H2 " 
		elif(line_edit[12:16]==" H22"): new_atom_name="2H2 " 
		elif(line_edit[12:16]==" H41"): new_atom_name="1H4 " 
		elif(line_edit[12:16]==" H42"): new_atom_name="2H4 " 
		elif(line_edit[12:16]==" H61"): new_atom_name="1H6 " 
		elif(line_edit[12:16]==" H62"): new_atom_name="2H6 " 

		elif(line_edit[12:16]==" OP1"): new_atom_name=" O1P"
		elif(line_edit[12:16]==" OP2"): new_atom_name=" O2P"

		if(	new_atom_name!=""): line_edit=line_edit[0:12]+new_atom_name+line_edit[16:]

		if(line_edit.count('OP1')> 0): error_exit_with_message("line_edit.count('OP1')> 0")
		if(line_edit.count('OP2')> 0): error_exit_with_message("line_edit.count('OP2')> 0")


		###################################################################################################################

		outid.write(line_edit)


fastaid.write('\n')
outid.close()
  #fastaid.close()

