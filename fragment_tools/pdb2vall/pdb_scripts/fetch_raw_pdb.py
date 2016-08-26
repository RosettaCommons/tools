#!/usr/bin/env python
from sys import argv
from os import system,path
import os
import ConfigParser

labdatabases = os.environ["PDB_DIR"]

# remote host for downloading pdbs
remote_host = os.environ["INET_HOST"]

pdbname = argv[1]

netpdbname = labdatabases
netpdbhost = ""
if remote_host:
    netpdbhost = remote_host

if not netpdbname.endswith(path.sep):
    netpdbname += path.sep

netpdbname += pdbname[1:3] + '/' + pdbname

if netpdbhost:
    system("scp %s:%s.pdb ." % ( netpdbhost, netpdbname) )
else:
    system("cp %s.pdb ." % netpdbname )
