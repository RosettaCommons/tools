#!/bin/bash


if [ -z "$1" ]
then
    echo Usage: $(basename $0) pdb_code_or_ALL \[rms=5\]
    echo examine the dwindling of close \(\<rms\) structures through the filters
    exit
fi
   
if [ -z "$2" ]
then
    rms=5
else
    rms=$2
fi
    
basedir=$PWD

if [ "$1" = "ALL" ]
then
    pdbs=$('ls' -d ????)
else
    pdbs=$1
fi

printf "Number of decoys within %i Angstroms in " $rms
echo $pdbs
printf "%10s %10s %10s %10s %10s %10s\n" "pdb" "basis" "all" "decoys" "top5k" "debumped"

for pdb in $pdbs
do

    scorefiles=$basedir/$pdb/scorefiles
    if [ ! -x $scorefiles ]; then scorefiles=$basedir/scorefiles; fi
    if [ ! -x $scorefiles ]
    then
	printf "%10s scorefiles not found\n" $pdb
    else

	cd $scorefiles

	printf "%10s" $pdb
	basis=$(grep -vE 'filename|nat|inp' $pdb.sc | wc | gawk {'print $1;'})
	printf " %10d" $basis
	for ext in sc decoys.sc top5000.sc top5000.nobumps.sc 
	do 
	    n=$(grep -vE 'nat|inp' $pdb.$ext | \
		filter_on.pl rms lt $rms 2>/dev/null | \
		grep -v 'filename' | wc | gawk {'print $1;'})
	    printf " %10d" $n
	done
	
	printf "\n"
	cd $basedir

    fi

done