#!/bin/bash

# summarize the clusters of each pdb in the current directory
# JJG 12/5/1

pdbs=$('ls' -d ????)
setname=$(basename $PWD)

echo Set $setname
date
echo --------------------------------------------

for pdb in $pdbs
do
    scorefiles=$pdb/scorefiles
    echo $pdb:
    echo
    head $pdb/top200/clusterscores.bysize
    printf "         %19s %7.2f\n" "native" \
	$(grep native $scorefiles/$pdb.native.sc | gawk {'print $2'})
    printf "         %19s %7.2f %6.2f\n" "nat_mcm" \
	$(grep nat_mcm $scorefiles/$pdb.native.sc | gawk {'print $2 " " $3'})
    printf "         %19s %7.2f %6.2f\n" "inp_rep" \
	$(grep inp_rep $scorefiles/$pdb.native.sc | gawk {'print $2 " " $3'})
    echo --------------------------------------------
done
