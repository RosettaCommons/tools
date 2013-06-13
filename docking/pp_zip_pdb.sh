#!/bin/bash


if [ -z "$1" ]
then
    echo Usage: $(basename $0) pdb_code_or_ALL
    echo zip the two-letter directories of pdb files
    exit
fi
    
basedir=$PWD

if [ "$1" = "ALL" ]
then
    pdbs=$('ls' -d ????)
    echo zipping ALL: $pdbs
else
    pdbs=$1
fi

for pdb in $pdbs
do

    cd $basedir/$pdb
    
    subdirs=$('ls' -d ??)
    nsubdirs=$(echo $subdirs | wc | gawk '{print $2}')
    echo zipping the $nsubdirs directories of structures of $pdb
    
    for dir in $subdirs
    do
	echo zipping $pdb/$dir
	tar czf $dir.tgz $dir --remove-files
	rmdir $dir
    done
    
    cd $basedir

done