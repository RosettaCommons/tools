#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
import os
import copy

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options

##############################USAGE#################################################################

#USAGE:THE CODE WILL GENERATE TORSION OF THE WHOLE STRUCTURE. If sample_res is specified, anchor
#additional file contain the torsions of nts corresponding to the sample is also generated.

#extract_torsions_from_pdb.py -pdb 1zih_RNA.pdb -whole_structure True

#extract_torsions_from_pdb.py -pdb Rosetta_1jj2_NO_5S.pdb  -whole_structure True -skip_idealize True 

#extract_torsions_from_pdb.py -pdb expand_radius_50_1nuj_RNA_A.pdb -sample_res 8 9 10 11

#extract_torsions_from_pdb.py -pdb 1l2x_RNA_B_with_sym_partner.pdb  -sample_res 18 19 20 21 22 23 24

#extract_torsions_from_pdb.py  -pdb 1_30_1jj2_RNA.pdb -sample_res 5 6 7

#bsub -W 140:0 -M 4000000 -o log.out -e log.err -J extract_torsion_from_pdb_1NUJ extract_torsions_from_pdb.py -pdb 1jj2_RNA.pdb -sample_res 31 32 33


####################################################################################################


copy_argv=copy.deepcopy(argv)

#####################################
print "Enter %s " %(list_to_string(copy_argv) )

HOMEDIR = expanduser('~') 

print "HOMEDIR= %s" %(HOMEDIR)


#exe_folder=HOMEDIR + "/src/mini/bin/" 

exe_folder=HOMEDIR + "/src/NO_GRAPHIC_mini/bin/" #Switch to the no_graphic version on Nov 21, 2011

if( exists(exe_folder )==False):  error_exit_with_message("Cannot find exe_folder: %s " %(exe_folder) )

PARIN_TEST= exe_folder + "/parin_test.macosgccrelease" 
if( exists(PARIN_TEST)==False): PARIN_TEST = exe_folder + "/parin_test.linuxgccrelease"


RNA_DATABASE= exe_folder + "/rna_database.macosgccrelease" 
if( exists(RNA_DATABASE)==False): RNA_DATABASE = exe_folder + "/rna_database.linuxgccrelease"


database_folder= HOMEDIR + "/minirosetta_database/" 

if( exists(PARIN_TEST)==False):  error_exit_with_message("PARIN_TEST_EXE (%s) doesn't exist! " %(PARIN_TEST) )
if( exists(RNA_DATABASE)==False):  error_exit_with_message("RNA_DATABASE_EXE (%s) doesn't exist! " %(RNA_DATABASE) )
if( exists(database_folder )==False):  error_exit_with_message("database_folder (%s) doesn't exist!"   %(database_folder) )

#####################################################################



def extract_sample_res_torsions(sample_res_list, local_pdb_file):

	print "extract_sample_res_torsion for pdb=%s " %(local_pdb_file)
	######################################################################
	# rna_database.macosgccrelease  -database ~/minirosetta_database -vall_torsions -s idealize_2gxb_RNA.pdb  -o idealize_2gxb_RNA.torsions

	full_length_torsion_outfile="full_length_" + local_pdb_file.replace(".pdb", "") + ".torsions"
	sample_res_torsion_outfile="sample_res_" + local_pdb_file.replace(".pdb", "") + ".torsions"


	if(exists(full_length_torsion_outfile)): submit_subprocess("rm %s " %(full_length_torsion_outfile) )
	if(exists(sample_res_torsion_outfile)): submit_subprocess("rm %s " %(sample_res_torsion_outfile) )

	extract_pdb_torsion_command=RNA_DATABASE
	extract_pdb_torsion_command += " -vall_torsions " #/algorithm
	extract_pdb_torsion_command += " -database %s" %(database_folder)
	extract_pdb_torsion_command += " -s %s " %(local_pdb_file)
	extract_pdb_torsion_command += " -o %s " %(full_length_torsion_outfile) #output

	extract_pdb_torsion_command += " > extract_%s_torsion_log.out 2> extract_%s_torsion_log.err" %(local_pdb_file, local_pdb_file)
	print "extract_pdb_torsion_command=%s " %(extract_pdb_torsion_command)
	submit_subprocess( extract_pdb_torsion_command )


	#########################################################################
	#Excise out the portion of the pdb that corresponds to the native loop.

	if(len(sample_res_list)==0):

		return full_length_torsion_outfile

	else:

		SAMPLE_RES_TORSION = open( sample_res_torsion_outfile, 'w')

		full_length_torsion_list=open(full_length_torsion_outfile).readlines()

		print "extract_torsion_sample_res= ", sample_res_list
		print "(seq_num, include_line):"

		for n in range(len(full_length_torsion_list)):
			line=full_length_torsion_list[n]
			line_split=line.split()

			seq_num=int(line_split[-1])
			if(seq_num!=(n+1)): error_exit_with_message(	"seq_num=(%d)!=(%d)=(n+1)" %(seq_num, n+1)) #This ensures that PDB is renumbered.

			include_line=False

			if(seq_num in sample_res_list): include_line=True
			if( (seq_num-1) in sample_res_list): include_line=True
			if( (seq_num+1) in sample_res_list): include_line=True

			print "(%3d," %(seq_num) ,
			if(include_line):
				print "T",
			else:
				print "F",
			print ")",
		
			if(include_line): SAMPLE_RES_TORSION.write( line )

		print 
		SAMPLE_RES_TORSION.close()

		return sample_res_torsion_outfile


#########################################################################


skip_idealize=parse_options( argv, "skip_idealize", "False") #Nov 21, 2011

pdb_file=parse_options( argv, "pdb", "")

pdb_file=os.path.abspath(pdb_file)

if(exists(pdb_file)==False): error_exit_with_message("pdb_file (%s) doesn't exist!" %(pdb_file))

sample_res_list= parse_options( argv, "sample_res", [-1] )

if(len(sample_res_list)==0): 

	whole_structure= parse_options( argv, "whole_structure", "False" )

	if(whole_structure):
		print "WARNING, user did not pass in sample_res_list, calculate torsions for whole structure!"
	else:
		error_exit_with_message("len(sample_res_list)==0 and whole_structure==False")

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )

#####################################################################
main_folder=os.path.abspath(".")

foldername="extract_torsions_%s" %( basename(pdb_file).replace(".pdb","") )

if( exists(foldername) ):  
	print "foldername (%s) already exist! ...removing...." %(foldername)
	submit_subprocess("rm -r %s " %(foldername) )

submit_subprocess( 'mkdir %s' %(foldername)  )

os.chdir( foldername)

if(pdb_file==os.path.abspath(basename(pdb_file))): error_exit_with_message("pdb_file==os.path.abspath(basename(pdb_file))!")

submit_subprocess("cp -r %s %s " %(pdb_file, basename(pdb_file) ) )
 
pdb_file=basename(pdb_file)

submit_subprocess("renumber_pdb_in_place.py %s " %( pdb_file ) ) #important!

#######################create_idealize_RNA###################################
# rna_test.macosgccrelease  -database ~/minirosetta_database -s 2gxb_RNA.pdb  -rna_idealize

if(skip_idealize==False):
	idealize_pdb_file="idealize_%s" %(pdb_file)


	IDEALIZE_RNA_command=PARIN_TEST
	IDEALIZE_RNA_command += " -algorithm rna_idealize " #/algorithm
	IDEALIZE_RNA_command += " -database %s" %(database_folder)
	IDEALIZE_RNA_command += " -s %s " %(pdb_file)
	IDEALIZE_RNA_command += " > idealize_rna_log.out 2> idealize_rna_log.err"

	print "IDEALIZE_RNA= %s " %(IDEALIZE_RNA_command)
	submit_subprocess( IDEALIZE_RNA_command )

	if(exists(idealize_pdb_file)==False): error_exit_with_message("idealize_pdb_file (%s) doesn't exist!" %(idealize_pdb_file) )

if(skip_idealize==False):
	IDEALIZED_sample_res_torsion_outfile=extract_sample_res_torsions(sample_res_list, idealize_pdb_file) #Before Nov 21 2011, there was error where the last argument here is pdb_file | should be effect final result


NOT_IDEALIZED_sample_res_torsion_outfile=extract_sample_res_torsions(sample_res_list, pdb_file) #Before Nov 21 2011, there was error where the last argument here is idealize_pdb_file| should be effect final result


os.chdir( main_folder)

#############MOVE TO MAIN_FOLDER############################################
if(skip_idealize==False):
	if(exists(IDEALIZED_sample_res_torsion_outfile)): submit_subprocess("rm  %s " %(IDEALIZED_sample_res_torsion_outfile) )

if(exists(NOT_IDEALIZED_sample_res_torsion_outfile)): submit_subprocess("rm %s " %(NOT_IDEALIZED_sample_res_torsion_outfile) )

#############MOVE TO MAIN_FOLDER############################################

if(skip_idealize==False):
	submit_subprocess("mv %s/%s %s " %(foldername,IDEALIZED_sample_res_torsion_outfile, IDEALIZED_sample_res_torsion_outfile) )

submit_subprocess("mv %s/%s %s " %(foldername,NOT_IDEALIZED_sample_res_torsion_outfile, NOT_IDEALIZED_sample_res_torsion_outfile) )


print "Exit %s " %(list_to_string(copy_argv) )


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
