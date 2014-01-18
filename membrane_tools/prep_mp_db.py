#!/usr/bin/python
## @file prep_mpdb.py
##
## @brief 	Prep Data files for individual membrane proteins
## @details For a given membrane protein run, grabs tr pdb, splits by chain, 
##			generates a span file from structure, generates a lipid accessibility file, 
##			outputs a fasta file and generates and xml file based on the setup of these resources for
##			Rosetta membrane
##
## @author Rebecca Alford
## @note Last Modified: 1/17/14