#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
from glob import glob
import os
import copy

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
from SWA_dagman_python.utility.PDB_operations import *
##############################USAGE#################################################################

#create_torsion_database.py -pdb_folder ../../RNA09_WORKINGSRC/RNA09PDBs/ > CREATE_DB_LOG.out 2> CREATE_DB_LOG.err

#create_torsion_database.py -pdb_folder ../../RNA09_WORKINGSRC/SUBSET_TEST/

#extract_torsions_from_pdb.py -pdb 1zih_RNA.pdb -whole_structure True


####################################################################################################


copy_argv=copy.deepcopy(argv)

#####################################
print "Enter %s " %(list_to_string(copy_argv) )

HOMEDIR = expanduser('~') 

print "HOMEDIR= %s" %(HOMEDIR)

#############################################

pdb_folder=parse_options( argv, "pdb_folder", "")

if(exists(pdb_folder)==False): error_exit_with_message("pdb_folder (%s) doesn't exist!")

pdb_folder=os.path.abspath(pdb_folder)

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )

#####################################################################
main_folder=os.path.abspath(".")

foldername="create_torsion_DB_%s" %(basename(pdb_folder))

if( exists(foldername) ):  
	print "foldername (%s) already exist! ...removing...." %(foldername)
	submit_subprocess("rm -r %s " %(foldername) )

submit_subprocess( 'mkdir %s' %(foldername)  )

os.chdir( foldername)

globstring = pdb_folder+'/*.pdb'
print "globstring= ", globstring

globfiles = glob( globstring )
globfiles.sort()
if(len( globfiles ) == 0): error_exit_with_message("len( globfiles ) == 0!")

PDB_src_table="PDB_src_table.txt"

PDB_SRC_TABLE=open(PDB_src_table, 'w')


PDB_SRC_COLUMN_LIST=["PDB            ", "LENGTH    ", "TITLE     "]

for n in range(len(PDB_SRC_COLUMN_LIST)):
	PDB_SRC_TABLE.write("%-*s" %(len(PDB_SRC_COLUMN_LIST[n])+10, PDB_SRC_COLUMN_LIST[n]) )

PDB_SRC_TABLE.write('\n')

num_RNA_struct=0
total_nts=0

all_torsions_file="RICHARDSON_RNA09.torsions"

if(exists(all_torsions_file)): error_exit_with_message("all_torsions_file (%s) already exist!")

for pdb_file in globfiles:

	if( os.path.abspath(pdb_file)==os.path.abspath(basename(pdb_file)) ): error_exit_with_message("os.path.abspath(pdb_file)=os.path.abspath(basename(pdb_file))")

	submit_subprocess("cp -r %s %s " %(pdb_file, basename(pdb_file) ) )
 
	pdb_file=basename(pdb_file)

	######get the title of the PDB file######
	PDB_LINES = open( pdb_file ).readlines()
	title=""
	for pdb_line in PDB_LINES:
		if(pdb_line[0:5]=="TITLE"): title+=pdb_line[10:-1]

	if(len(title)>100): title=title[:100] + " CONT."

	'''
	TITLE     THREE-DIMENSIONAL STRUCTURE OF RIBONULCEASE T1 COMPLEXED
	TITLE    2 WITH AN ISOSTERIC PHOSPHONATE ANALOGUE OF GPU: ALTERNATE
	TITLE    3 SUBSTRATE BINDING MODES AND CATALYSIS.
	01234567890

	CRYSTAL STRUCTURE OF SELENIUM-MODIFIED DIELS-ALDER RIBOZYME COMPLEXED WITH THE PRODUCT OF THE REACTION BETWEEN N- PENTYLMALEIMIDE AND COVALENTLY ATTACHED 9-
	0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
	0         1         2         3         4         5         6         7         8         9         10         11       12        13        14        15
	'''
	#########################################

	start_pdb_name=pdb_file

	if(exists(pdb_file)==False): error_exit_with_message("pdb_file (%s) doesn't exist!" %(pdb_file))

	make_rosetta_ready_command="SWA_make_rna_rosetta_ready.py %s -output_pdb %s " %(pdb_file, "rs_ready_" + pdb_file)

	submit_subprocess(make_rosetta_ready_command)

	pdb_file="rs_ready_" + pdb_file

	if(exists(pdb_file)==False): error_exit_with_message("pdb_file (%s) doesn't exist!" %(pdb_file))

	submit_subprocess("Rosetta_import_and_dump_pdb.py -s %s" %(pdb_file))

	pdb_file="rosetta_%s" %(pdb_file)

	if(exists(pdb_file)==False): error_exit_with_message("pdb_file (%s) doesn't exist!" %(pdb_file))

	replace_chain_inplace_func(pdb_file, "A") #submit_subprocess("replace_chain_inplace.py %s A" %(current_pdb) )
	renumber_pdb_in_place_func(pdb_file)			 #submit_subprocess("renumber_pdb_in_place.py %s" %(current_pdb) )


	######Extract the torsions#######################

	extract_torsion_command="extract_torsions_from_pdb.py -pdb %s -whole_structure True -skip_idealize True" %(pdb_file)

	submit_subprocess(extract_torsion_command)

	full_length_torsion_outfile="full_length_" + pdb_file.replace(".pdb", "") + ".torsions"

	if(exists(full_length_torsion_outfile)==False): error_exit_with_message("full_length_torsion_outfile (%s) doesn't exist!")

	submit_subprocess("cat %s >> %s " %(full_length_torsion_outfile, all_torsions_file))

	######Get the number of nts in the pdb_file ######
	
	fasta_file=os.path.abspath("./fasta_%s" %(pdb_file))
	submit_subprocess('SWA_pdb2fasta.py %s > %s' %(pdb_file,fasta_file))

	SEQUENCE = (open( fasta_file ).readlines()[1][:-1]).upper() 
	RNA_length=len(SEQUENCE)

	############Check that the nts length match that in the torsion_file######
	torsion_file_length=len(open( full_length_torsion_outfile ).readlines())

	if(RNA_length!=torsion_file_length): error_exit_with_message("RNA_length=(%s)!=(%s)=torsion_file_length for pdb_file (%s)" %(RNA_length, torsion_file_length, pdb_file))

	##########################################################################

	total_nts+=RNA_length
	num_RNA_struct+=1

	PDB_SRC_TABLE.write("%-*s" %(len(PDB_SRC_COLUMN_LIST[0])+10, start_pdb_name) )
	PDB_SRC_TABLE.write("%-*s" %(len(PDB_SRC_COLUMN_LIST[1])+10, RNA_length) )
	PDB_SRC_TABLE.write("%-*s" %(len(PDB_SRC_COLUMN_LIST[2])+10, title) )


	PDB_SRC_TABLE.write('\n')
	PDB_SRC_TABLE.flush()


summary_line= "SUMMARY: NUM_RNA_STRUCT=%s | TOTAL_NTS=%s" %(num_RNA_struct, total_nts)

print summary_line
PDB_SRC_TABLE.write(summary_line + '\n') 


PDB_SRC_TABLE.close()


print "-------------------------------------------------------------------------------------------"
print "Successfully RAN: %s" %(list_to_string(copy_argv))
print "-------------------------------------------------------------------------------------------"



'''
####################EMAIL TO TONY################################

Hi Tony!

Edit: 

src/pilot_apps.src.settings.all

or

src/apps.src.settings

These files list all the main executables to compile. Make sure at leas t one of them include:

 "pilot/rhiju" : [ "rna_test" ]

and

 "public/rna" : [ "rna_database" ]


Cheers,
Rhiju

Examples below (subsitute 2gxb with your favorite pdb!):

First, idealize it:

 rna_test.macosgccrelease  -database ~/minirosetta_database -s 2gxb_RNA.pdb  -rna_idealize

Then get the torsions:

 rna_database.macosgccrelease  -database ~/minirosetta_database -vall_torsions -s idealize_2gxb_RNA.pdb  -o idealize_2gxb_RNA.torsions

When you run fragment assembly, then include the command-line '-vall_torsions idealize_2gxb_RNA.torsions'


        Cheers,
         Rhiju


####################EMAIL TO PARIN################################

For native torsions, it is:

rna_database.<exe> -vall_torsions -s <input pdb> -o <outfile> -database <database>

To idaelize ahead of time, I think you can use:

rna_test.<exe> -idealize   [or maybe is it is -idealize_rna or something, sorry]

I typically have extracted torsions from both the starting pdb and the idealized starting pdb, and concatenated the two files.

Cheers,
 Rhiju

##############################################################################


rna_test.<exe> -rna_idealize 
'''
