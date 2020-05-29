#!/usr/bin/python

import numpy as np
import argparse
import os.path
from recces.util import compute_rigid_body_ref

parser = argparse.ArgumentParser( description = "compute reference energy for base pair energetics" )
parser.add_argument('RMSD_cutoff',type=float,nargs='?',help='RMSD cutoff for 6D rot/translations (Angstroms)',default=2.0 )
parser.add_argument('xyz_file',nargs='?',help='txt file with xyz atoms for object, centered at origin',default='xyz.txt')


args = parser.parse_args()

if not os.path.isfile( args.xyz_file ):
    print "Please specify a valid file with -xyz_file. Could not find ", args.xyz_file
    exit( 0 )

ref_energy = compute_rigid_body_ref( args.RMSD_cutoff, args.xyz_file )

print 'Analytical reference energy at 1 M for RMSD cutoff %4.1f with xyz from %s: %f' % ( args.RMSD_cutoff, args.xyz_file, ref_energy )

