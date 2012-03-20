#!/usr/local/bin/python
from sys import argv
from os import system,path
import ConfigParser

## get PDB2VALL_PATH
PDB2VALL_PATH = path.abspath(path.dirname(__file__))
PDB2VALL_PATH = PDB2VALL_PATH[0:PDB2VALL_PATH.rfind("/pdb2vall/")];
PDB2VALL_PATH += "/pdb2vall/"
if not PDB2VALL_PATH:
    stderr.write("ERROR: you should specify the path where your packages are first.\n")
    exit()

## read config
config = ConfigParser.RawConfigParser()
config.read(PDB2VALL_PATH + "pdb2vall.cfg")

pdbname = argv[1]

netpdbname = config.get('pdb2vall', 'wwpdb_path')
netpdbhost = ""
if config.has_option('pdb2vall', 'wwpdb_host'):
    netpdbhost = config.get('pdb2vall', 'wwpdb_host')

if not netpdbname.endswith(path.sep):
    netpdbname += path.sep

netpdbname += pdbname[1:3] + '/' + pdbname

if netpdbhost:
    system("scp %s:%s.pdb ." % ( netpdbhost, netpdbname) )
else:
    system("cp %s.pdb ." % netpdbname )
