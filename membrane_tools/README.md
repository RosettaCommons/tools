Rosetta Membrane Tools
====================
This directory contains scripts for preping and running Rosetta membrane applications with the mew membrane protien framework. Below are some quick descriptions of each script - more detailed information can be found on the Rosetta docuemntation Wiki: 
* run_lips.pl - Generate a lipophobicity data file given a .span and .fasta file
* alignblast.pl - Dependency required for run_lips.pl
* mptest_ut.py - Run all unit tests associated with the membrane framework
* octopus2span.pl - Convert the .oct data file from OCTOPUS to a .span file format
* prep_mp_db.py - Prep required data files for a given PDB ID (only transformed pbds please!)   
* write_mp_xml.py - Write mmebrane resource definition file for running the framework