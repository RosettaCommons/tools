#!/bin/bash
#Author: SAQ
#email: ssaq@jhu.edu

if [ ! -e "$1" ]; then
        echo "Usage: ./mini_dock.bash <pdb>"
	echo "Do not add .pdb extension to the pdb ID"
        exit
fi

#input pdb ID
pdb=$1

#unique folder created to store the decoys. Change this if you have more than one run in the same directory
#must be two letters
prefix=aa

#number of decoys (test = 20, standard = 1000)
nstruct=1000



#--------
#standard local perturbation run
protocol_flags="-dock_pert 3 8 -spin -dock_mcm"

#standard global run
#protocol_flags="-randomize1 -randomize2 -spin -dock_mcm"
#--------


#if you want to fix the sidechains of partner 1, add -norepack1
#if you want to fix the sidechains of partner 2, add -norepack2
sidechain_flags="-ex1 -ex2aro"

mkdir $pdb
mkdir $pdb/$prefix

$minirosetta_bin/docking_protocol.linuxgccrelease -database $minirosetta_database -out:prefix $prefix -s ./prepack/$pdb.ppk.pdb $protocol_flags $sidechain_flags -nstruct $nstruct -mute core -no_optH -out:file:fullatom -out:path:pdb $pdb/$prefix -out:file:o $pdb

mv $prefix"$pdb".fasc ./$pdb 
