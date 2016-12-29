#!/usr/bin/env python
'''
In this example, pdb files are converted to Amber coordinate files (rst7). One topology file (parm7) is created for both rst7 files, as the PDBs have the same amino acid sequence.
The rst7 files are minimized using Amber's sander minimization engine, and the energies of the minimized structures are gathered and printed.
'''

import AmbRose
from glob import glob
import os
import pytraj as pt

os.chdir('example_input')
pdbfiles = glob("*.pdb")

## Convert the pdbs to Amber rst7 files and a parmfile.
## Arguments: [0] = list of pdbfiles,
##            [1] = make one parm file for all rst7 files? ('0' will make a parmfile for each rst7, '1' will make one parmfile called 'parm.parm7')
rst7files, parmfiles = AmbRose.convert_pdbs_to_rst7_parm7_files( pdbfiles, 1 )

#print(rst7files)
#print(parmfiles)

## Minimized the rst7 files under one parm file.
## Arguments: [0] = parmfile (str), 
##            [1] = list of rst7files, 
##            [2] = number of cores to use for MPI minimization (int), 
##            [3] = overwrite? (0 or 1)
AmbRose.batch_minimize_with_sander( parmfiles[0], rst7files, 1, 0 )

## Gather the minimized rst7 files into a list.
minimized_rst7s = glob("min*.rst7")

print "Making trajectory"
minimized_traj = pt.iterload(minimized_rst7s, parmfiles[0])

print "Getting Amber Energies"
## get_energies takes the trajectory and the list of rst file names
energy_data = AmbRose.get_energies( minimized_traj, minimized_rst7s )

AmbRose.print_energies( energy_data )

## write_energies takes energy data and an outfile name.
AmbRose.write_energies( energy_data, 'scores.sc' )
