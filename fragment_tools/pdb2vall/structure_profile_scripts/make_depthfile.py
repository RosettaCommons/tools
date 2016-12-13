#!/usr/bin/env python

from sys import argv, exit, stderr
from os import system, getcwd, chdir, path
from os.path import exists
import ConfigParser

## 111206 fixed a bug where the picker doesn't read the results from multiple chains -
##        modify that into running program with multiple chains, and parse the single chain results out.

## GET PDB AND RUN DEPTH
if len( argv ) < 2:
    print
    print 'USAGE: %s 2oxgA.pdb' % argv[0]
    print
    exit()

pdbname  = argv[1][:4]
pdbchain = argv[1][4]

## get PDB2VALL_PATH
PDB2VALL_PATH = path.abspath(path.dirname(__file__))
PDB2VALL_PATH = PDB2VALL_PATH[0:PDB2VALL_PATH.rfind("/pdb2vall/")] + "/";
if not PDB2VALL_PATH:
    stderr.write("ERROR: should specify the path of the pdb2vall package.\n"); exit()

## read config
config = ConfigParser.RawConfigParser(allow_no_value=True)
config.read(PDB2VALL_PATH + "pdb2vall/pdb2vall.cfg")

# default
DEPTH = PDB2VALL_PATH + "pdb2vall/structure_profile_scripts/DEPTH-CLONE-2.8.7/DEPTH";
if config.get('pdb2vall', 'depth'):
    DEPTH = config.get('pdb2vall', 'depth')

DEPTH_THREADS = config.get('pdb2vall', 'depth_num_cpus')
if not DEPTH_THREADS:
    DEPTH_THREADS = 1

## MAKE A SUBDIR TO RUN THAT
basedir = getcwd() + "/"
system("mkdir %s" %( pdbname+pdbchain ))
chdir( basedir + pdbname+pdbchain )
system("rm *.log");

## GET WHOLE CHAIN PDB AND CLEAN THAT UP.
cmd = PDB2VALL_PATH + 'pdb2vall/pdb_scripts/fetch_raw_pdb.py %s' % pdbname
system( cmd )

print 'cleaning up pdb'
system( PDB2VALL_PATH + 'pdb2vall/pdb_scripts/clean_pdb.py %s.pdb > /dev/null' % pdbname )
system('rm -f %s.pdb' % argv[1][:4] )
system('ln -s %s.pdb.clean.pdb %s.pdb' % ( pdbname, pdbname ))

print 'running DEPTH'
cmd_depth = DEPTH + ' -thread ' + DEPTH_THREADS + ' -i %s.pdb -o %s.depth' % ( pdbname, pdbname )
system( cmd_depth )

wholechain_results  = pdbname + '.depth-residue.depth'
singlechain_results = pdbname+pdbchain + '.depth-residue.depth'

if exists( wholechain_results ):
    cmd = "head -1 %s > %s ; grep ^%s %s >> %s" %( wholechain_results, singlechain_results, pdbchain, wholechain_results, singlechain_results )
    system( cmd )

if not exists( wholechain_results ):
    print 'fail to get the depth data from running the whole chains pdb - now try running with pdb with only a chain'

    cmd = PDB2VALL_PATH + 'pdb2vall/pdb_scripts/get_pdb_new.py %s %s' % ( pdbname, pdbchain )
    system( cmd )

    print 'running DEPTH'
    cmd_depth = DEPTH + ' -thread ' + DEPTH_THREADS + ' -i %s.pdb -o %s.depth' % ( pdbname+pdbchain, pdbname+pdbchain )
    system( cmd_depth )

if exists( singlechain_results ):
    chdir( basedir )
    system("ln -s ./%s/%s ." %( pdbname+pdbchain, singlechain_results ))
else:
    print "job failed - something bad happened here"

