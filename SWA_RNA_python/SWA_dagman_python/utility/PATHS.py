#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser,abspath
from sys import exit, argv
import string
from time import sleep
import glob
import time
import fileinput
import os

######################################################################
from SWA_dagman_python.utility import USER_PATHS

from SWA_dagman_python.utility.error_util import *
from SWA_dagman_python.utility.list_operations import *

##################NOTES###############################################
#Current the scripts should be compatible with either a MacOS or Linux system
#Using python version 2.3 or above
#NOT yet tested and probably won't work with Windows. For example currently uses SLASH (/) for directory seperation and not BACK-SLASH (for Windows)

######################################################################
HOMEDIR = expanduser('~')


RELEASE_MODE=True
######################################################################
if(RELEASE_MODE==True):
	ROSETTA_BIN_FOLDER             = os.path.expanduser(USER_PATHS.USER_ROSETTA_BIN_FOLDER)
	NO_GRAPHIC_ROSETTA_BIN_FOLDER  = os.path.expanduser(USER_PATHS.USER_ROSETTA_BIN_FOLDER)
	ROSETTA_DATABASE_FOLDER        = os.path.expanduser(USER_PATHS.USER_ROSETTA_DATABASE_FOLDER)
	NEW_SRC_CODE                   = True
	NO_GRAPHIC                     = True

else:



	ROSETTA_BIN_FOLDER             = HOMEDIR + "/ROSETTA/rosetta_CURRENT/rosetta_source/bin/"
	NO_GRAPHIC_ROSETTA_BIN_FOLDER  = HOMEDIR + "/ROSETTA/rosetta_CURRENT/NO_GRAPHIC_rosetta_source/bin/"
	ROSETTA_DATABASE_FOLDER        = HOMEDIR + "/ROSETTA/rosetta_CURRENT/rosetta_database/"
	NEW_SRC_CODE                   = False
	NO_GRAPHIC                     = False

	'''

	ROSETTA_BIN_FOLDER             = HOMEDIR + "/ROSETTA/rosetta_UPDATE_VERSION_2/rosetta_source/bin/"
	NO_GRAPHIC_ROSETTA_BIN_FOLDER  = HOMEDIR + "/ROSETTA/rosetta_UPDATE_VERSION_2/rosetta_source/bin/"
	ROSETTA_DATABASE_FOLDER        = HOMEDIR + "/ROSETTA/rosetta_UPDATE_VERSION_2/rosetta_database/"
	NEW_SRC_CODE                   = True
	NO_GRAPHIC                     = False

	'''

######################################################################

#This assumes that as PATHS.py in located in a folder one level deeper than the main python folder.
SWA_DAGMAN_PYTHON_MAIN_FOLDER=list_to_string(__file__.split('/')[:-2],'/')[1:]+'/'

SWA_dagman_python_PATHS=['', 'SWA_DAG/', 'FARFAR_DAG/', 'dagman/', 'misc/']
if(RELEASE_MODE==False): SWA_dagman_python_PATHS.extend(['chemical_shift/', 'projects/'])

SWA_dagman_python_PATHS=		map( lambda rel_path  : SWA_DAGMAN_PYTHON_MAIN_FOLDER + rel_path , SWA_dagman_python_PATHS )

Existing_PATHS=os.environ["PATH"].split(os.pathsep)


Do_update_PATHS=False

for n in range(len(SWA_dagman_python_PATHS)):

	if(n >= len(Existing_PATHS)):
		Do_update_PATHS=True
		break

	if(Existing_PATHS[n]!=SWA_dagman_python_PATHS[n]):
		Do_update_PATHS=True
		break

if(Do_update_PATHS):

	#print "BEFORE os.environ[\"PATH\"]=", os.environ["PATH"]

	os.environ["PATH"]=os.pathsep.join(SWA_dagman_python_PATHS)

	for n in range(len(Existing_PATHS)):

		if(Existing_PATHS[n] in SWA_dagman_python_PATHS): continue ##CHECK IF ALREADY EXIST!

		os.environ["PATH"]+=os.pathsep + Existing_PATHS[n]

	#print "AFTER os.environ[\"PATH\"]=", os.environ["PATH"]



##################################################################################################
#print "(PATHS.py) __file__.=         ", __file__
#print "SWA_DAGMAN_PYTHON_MAIN_FOLDER=", SWA_DAGMAN_PYTHON_MAIN_FOLDER
#print "ROSETTA_BIN_FOLDER=           ", ROSETTA_BIN_FOLDER
#print "NO_GRAPHIC_ROSETTA_BIN_FOLDER=", NO_GRAPHIC_ROSETTA_BIN_FOLDER
#print "ROSETTA_DATABASE_FOLDER=      ", ROSETTA_DATABASE_FOLDER
##################################################################################################

POSSIBLE_ROSETTA_EXES=[]

if(NEW_SRC_CODE):

	POSSIBLE_ROSETTA_EXES=['swa_rna_main', 'swa_rna_util', 'rna_extract', 'rna_denovo']
	if(RELEASE_MODE==False): POSSIBLE_ROSETTA_EXES.extend(['parin_test'])

else:
	POSSIBLE_ROSETTA_EXES=['rna_swa_test', 'rna_extract', 'rna_denovo']
	if(RELEASE_MODE==False): POSSIBLE_ROSETTA_EXES.extend(['parin_test'])

##################################################################################################

def PATH_exists(path_name):

	return exists(os.path.expanduser(path_name))

##################################################################################################
def use_new_src_code():

	return NEW_SRC_CODE

##################################################################################################
def is_release_mode():
	"""Release_mode specifies that features have been tested and
	is release for usage by the wider community."""

	return RELEASE_MODE

##################################################################################################

def SWA_unexpanduser(path_name_expanduser):
	#The reason why we need to unexpand the path is because while testing, I usually perform SWA and FARFAR setup locally and then rsync to cluster for job submission

	if(exists(path_name_expanduser)==False): error_exit_with_message("path_name_expanduser (%s) doesn't exists!" %(path_name_expanduser))

	if(path_name_expanduser[:len(HOMEDIR)]!=HOMEDIR):
		error_exit_with_message("path_name_expanduser[:len(HOMEDIR)]=(%s)!=(%s)=HOMEDIR" %(path_name_expanduser[:len(HOMEDIR)], HOMEDIR))

	path_name_with_tilde='~' + path_name_expanduser[len(HOMEDIR):]

	if(os.path.expanduser(path_name_with_tilde)!=path_name_expanduser): error_exit_with_message("os.path.expanduser(path_name_with_tilde)!=path_name_expanduser")

	if(PATH_exists(path_name_with_tilde)==False): error_exit_with_message("PATH_exist(path_name_with_tilde))==False")

	return path_name_with_tilde

##################################################################################################

def get_rosetta_EXE(EXE_name, no_graphic=NO_GRAPHIC):   #EXE=get_rosetta_EXE("parin_test", no_graphic)

	###PROBLEM IS THAT WHILE I PERSONALLY USE THE GRAPHIC MODE AS DEFAULT, RELEASE VERSION USUALLY USE THE NO_GRAPHIC!...SO:
	###Should make code flexible, in not if no_graphic=False, but if the graphic executable doesn't exist, then run no_graphic version (LOCATED IN THE SAME BIN FOLDER!)
	###REASON IS THAT some code like SWA_cluster.py will pass the no_graphic=True as default!###
	###ALSO SHOULD CHANGE SO THAT IF no no_graphic option is set to specificied into the script, then just call get_rosetta_EXE("exe_name") without the second variable!

	if(EXE_name not in POSSIBLE_ROSETTA_EXES): error_exit_with_message("Invalid EXE_name (%s) " %(EXE_name) )

	exe_folder=ROSETTA_BIN_FOLDER

	if(no_graphic): exe_folder=NO_GRAPHIC_ROSETTA_BIN_FOLDER

	#exe_folder=SWA_unexpanduser(exe_folder) # what?

	if( PATH_exists(exe_folder)==False):  error_exit_with_message("Cannot find exe_folder: %s " %(exe_folder) )

	LINUX_EXE = "%s%s.linuxgccrelease" %(exe_folder, EXE_name)
	LINUX_EXE2 = "%s%s.linuxclangrelease" %(exe_folder, EXE_name)

	MACOS_EXE = "%s%s.macosgccrelease" %(exe_folder, EXE_name)

	possible_EXE_count=0

	if(PATH_exists(LINUX_EXE)): possible_EXE_count+=1
	if(PATH_exists(LINUX_EXE2)): possible_EXE_count+=1

	if(PATH_exists(MACOS_EXE)): possible_EXE_count+=1

	if(possible_EXE_count!=1):

		print "possible_EXE_count=(%s)=!=1 for EXE_name (%s)" %(possible_EXE_count, EXE_name)
		##########################################################
		if(PATH_exists(LINUX_EXE)):
			print "LINUX_EXE (%s) exist!" %(LINUX_EXE)
		else:
			print "LINUX_EXE (%s) doesn't exist!" %(LINUX_EXE)

		if(PATH_exists(LINUX_EXE2)):
			print "LINUX_EXE2 (%s) exist!" %(LINUX_EXE2)
		else:
			print "LINUX_EXE2 (%s) doesn't exist!" %(LINUX_EXE2)

		if(PATH_exists(MACOS_EXE)):
			print "MACOS_EXE (%s) exist!" %(MACOS_EXE)
		else:
			print "MACOS_EXE (%s) doesn't exist!" %(MACOS_EXE)
		##########################################################

		error_exit_with_message("possible_EXE_count=(%s)=!=1 for EXE_name (%s)" %(possible_EXE_count, EXE_name))

	if(PATH_exists(LINUX_EXE)):  return LINUX_EXE
	if(PATH_exists(LINUX_EXE2)):  return LINUX_EXE2

	if(PATH_exists(MACOS_EXE)):  return MACOS_EXE

	error_exit_with_message("Should not reach this point of the function!")

######Since the exist condition checks below are alway called, will there slow down codes that call SWA_Util.py very often?######################
if(exists(ROSETTA_BIN_FOLDER)           ==False): error_exit_with_message("ROSETTA_BIN_FOLDER            (%s) doesn't exist" %(ROSETTA_BIN_FOLDER))
#if(exists(NO_GRAPHIC_ROSETTA_BIN_FOLDER)==False): error_exit_with_message("NO_GRAPHIC_ROSETTA_BIN_FOLDER (%s) doesn't exist" %(NO_GRAPHIC_ROSETTA_BIN_FOLDER))
if(exists(ROSETTA_DATABASE_FOLDER)      ==False): error_exit_with_message("ROSETTA_DATABASE_FOLDER       (%s) doesn't exist" %(ROSETTA_DATABASE_FOLDER))
if(exists(SWA_DAGMAN_PYTHON_MAIN_FOLDER)==False): error_exit_with_message("SWA_DAGMAN_PYTHON_MAIN_FOLDER (%s) doesn't exist" %(SWA_DAGMAN_PYTHON_MAIN_FOLDER))

#check that all the executable exist! OK comment this out, since for certain task, only some of the executables will be needed!
#for n in range(len(POSSIBLE_ROSETTA_EXES)):

#	if(ROSETTA_BIN_FOLDER!=""):
#		get_rosetta_EXE(POSSIBLE_ROSETTA_EXES[n], no_graphic=False)  #This checks that the EXE exist!

#	if(NO_GRAPHIC_ROSETTA_BIN_FOLDER !=""):
#		get_rosetta_EXE(POSSIBLE_ROSETTA_EXES[n], no_graphic=True) #This checks that the EXE exist!

##################################################################################################
def get_rosetta_EXE_specified_no_graphic_string(EXE_name, no_graphic_string):

	if(no_graphic_string==""): #no no_graphic_string specified
		return get_rosetta_EXE(EXE_name)
	else:
		if(no_graphic_string=="True"):
			print "Overwriting default no_graphic option to True, for EXE_name (%s)!" %(EXE_name)
			return get_rosetta_EXE(EXE_name, True)
		elif(no_graphic_string=="False"):
			print "Overwriting default no_graphic option to False, for EXE_name (%s)!" %(EXE_name)
			return get_rosetta_EXE(EXE_name, False)
		else:
			error_exit_with_message("Invalid no_graphic_test (%s)" %(no_graphic_string))

##################################################################################################

def get_rosetta_database_folder():

	if(exists( ROSETTA_DATABASE_FOLDER )==False ): error_exit_with_message("ROSETTA_DATABASE_FOLDER (%s) doesn't exist!" %(ROSETTA_DATABASE_FOLDER))

	#rosetta_database_with_tilde=SWA_unexpanduser(ROSETTA_DATABASE_FOLDER)

	#if(PATH_exists( rosetta_database_with_tilde )==False ): error_exit_with_message("rosetta_database_with_tilde (%s) doesn't exist!" %(rosetta_database_with_tilde))

	#return rosetta_database_with_tilde
	return ROSETTA_DATABASE_FOLDER

##################################################################################################

def get_PYDIR():

	#pydir_with_tilde=SWA_unexpanduser(SWA_DAGMAN_PYTHON_MAIN_FOLDER)

	#if(PATH_exists(pydir_with_tilde)==False): error_exit_with_message("pydir_with_tilde (%s) doesn't exist" %(pydir_with_tilde))

	#return pydir_with_tilde

	return SWA_DAGMAN_PYTHON_MAIN_FOLDER

##################################################################################################
def get_PYPATH(rel_path):

	if(isinstance( rel_path, str )==False):
		print "rel_path=", rel_path
		error_exit_with_message("rel_pathis not a instance of str!")

	#EXE_name could contain SLASH (/). For Windows, might want to convert SLASH (/) to BACK-SLASH

	pypath_with_tilde=get_PYDIR() + rel_path

	if(PATH_exists(pypath_with_tilde)==False): error_exit_with_message("pypath_with_tilde (%s) doesn't exist!" %(pypath_with_tilde))

	return pypath_with_tilde

##################################################################################################

def get_PYEXE(EXE_name):

	if(isinstance( EXE_name, str )==False):
		print "EXE_name=", EXE_name
		error_exit_with_message("EXE_name is not a instance of str!")

	#EXE_name could contain SLASH (/). For Windows, might want to convert SLASH (/) to BACK-SLASH

	if(get_PYDIR() in EXE_name or '~/' in EXE_name):
		pyexe_with_tilde=EXE_name
	else:
		pyexe_with_tilde=get_PYDIR() + EXE_name

	if(PATH_exists(pyexe_with_tilde)==False): error_exit_with_message("pyexe_with_tilde (%s) doesn't exist!" %(pyexe_with_tilde))

	return pyexe_with_tilde

##################################################################################################
'''
def get_PYDIR_with_tilde():

	pydir_with_tilde=SWA_unexpanduser(get_PYDIR())

	return pydir_with_tilde

##################################################################################################
def get_PYEXE_with_tilde(EXE_name):

	if(isinstance( EXE_name, str )==False):
		print "EXE_name=", EXE_name
		error_exit_with_message("EXE_name is not a instance of str!")

	#EXE_name could contain SLASH (/). For Windows, might want to convert SLASH (/) to BACK-SLASH

	pyexe_with_tilde=get_PYDIR_with_tilde()+EXE_name

	if(exists(os.path.expanduser(pyexe_with_tilde)==False)):
		error_exit_with_message("os.path.expanduser(pyexe_with_tilde) (%s) doesn't exist!" %(os.path.expanduser(pyexe_with_tilde)))

	return pyexe_with_tilde
'''



