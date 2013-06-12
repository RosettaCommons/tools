#!/bin/bash


if [ -z "$1" ]
then
    echo Usage: $(basename $0) pdb_code [rmscutoff=2.5]
    echo finishes post-processing a pdb that has been PUSHED to hep
    exit
fi

PATH=$PATH:~/simcode/scripts

pdb=$1  
rmscutoff=2.5
#if [ ! -z "$2" ]; then rmscutoff=5
if [ ! -z "$2" ]; then rmscutoff=$2
else 
	rmscutoff=5
	echo defaulting clustering rms cutoff to $rmscutoff; 
fi

pp_cluster_set.sh $pdb 10 $rmscutoff    | tee -a $pdb.pplog2.10
pp_cluster_set.sh $pdb 200 $rmscutoff   | tee -a $pdb.pplog2.200

if [ "$(basename $PWD)" == "unbound-E-7-16" ]; then
    pp_calc_contacts.bash $pdb | tee -a $pdb.pplog2.200
fi
