#!/bin/bash

# summarize the lowest scores of each pdb in the current directory
# JJG 7/9/2

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
    cutfa1=$(findIndex.pl trans $scorefiles/$pdb.top10.sc |gawk {'print $4'})
    cutfa2=$(findIndex.pl description $scorefiles/$pdb.top10.sc \
	     |gawk {'print $4'})
    head $scorefiles/$pdb.top10.sc | cut -c1-$cutfa1,$cutfa2-
    grep native $scorefiles/$pdb.native.sc | cut -c1-$cutfa1,$cutfa2-
    grep nat_mcm $scorefiles/$pdb.native.sc | cut -c1-$cutfa1,$cutfa2-
    grep inp_rep $scorefiles/$pdb.native.sc | cut -c1-$cutfa1,$cutfa2-
    echo --------------------------------------------
done
