#!/bin/bash
#Author: SAQ
#email: ssaq@jhu.edu

if [ ! -e "$1" ]; then
        echo "Usage: ./mini_prepack.bash <pdb>"
	echo "Do not add .pdb extension to the pdb ID"
        exit
fi




#input pdb ID 
pdb=$1

#Prepacking flags. used if docking is performed with rotamer trial minimization
protocol_flags="-dock_ppk"

#If you want to fix the sidechains of partner 1, add -norepack1
#If you want to fix the sidechains of partner 2, add -norepack2
sidechain_flags="-ex1 -ex2aro -use_input_sc"

mkdir prepack
cp $pdb.pdb ./prepack/

$minirosetta_bin/docking_protocol.linuxgccrelease -database $minirosetta_database -s $pdb.pdb $protocol_flags $sidechain_flags -nstruct 1 -mute core -no_optH -out:file:fullatom -out:path:pdb prepack -out:file:o $pdb

mv ./prepack/$pdb"_0001.pdb" ./prepack/$pdb.ppk.pdb 
rm $pdb.fasc 
