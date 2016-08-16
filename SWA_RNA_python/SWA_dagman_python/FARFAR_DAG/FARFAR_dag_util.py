#!/usr/bin/env python

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
import os
from os.path import basename, dirname, exists, expanduser
from time import sleep
from sets import Set
import copy
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################


def get_filtered_RMSD_file(DAG_ID):

	filename="DAG_ID_%d_filtered_RMSD.out" %(DAG_ID)
	return filename

def get_filtered_energy_file(DAG_ID):

	filename="DAG_ID_%d_filtered_energy.out" %(DAG_ID)
	return filename

def get_DAG_ID_folder(DAG_ID):

	foldername="DAG_ID_%d/" %(DAG_ID)
	return foldername

def get_SCORE_file(DAG_ID):

	filename="DAG_ID_%d_SCORE.out" %(DAG_ID)
	return filename

################################################

def get_FARFAR_ROSETTA_COMMAND_file():
	return "ROSETTA_COMMAND"

def get_FARFAR_common_args_file():
	return "FARFAR_common_args.txt"

def get_FARFAR_cluster_args_file():
	return "FARFAR_cluster_args.txt"

