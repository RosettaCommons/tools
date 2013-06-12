#!/bin/bash
#variable RMSD clustering

if [ -z "$3" ]
then
    echo Usage: $(basename $0) pdb_code n rmscutoff
    echo cluster a set of the top n structures of a given pdb with rmscutoff
    exit
fi

# if [ ! -x /usr/bin/R ] # commented out to use R in home directory
if [ ! -x /net/local/bin/R ]
then
    echo R not available on this workstation...no clustering possible
    echo consider pushing the files home with pp_push_set.sh
    exit
fi


basedir=$PWD
pdb=$1
ndecoys=$2
rmscutoff=$3
scorefiles=$basedir/$pdb/scorefiles

# assemble sets of decoys
echo clustering the top $ndecoys structures of $pdb in $basedir with a $rmscutoff A cutoff 

nplus1=$[ $ndecoys + 1 ]
topn=top$ndecoys

if [ ! -d "$scorefiles" ]; then echo No directory $scorefiles; exit; fi
cd $basedir/$pdb/$topn || exit

## do pairwise rms calculations
cat $scorefiles/$pdb.$topn | rms2.pl > rms.table

## calculate clusters from rms (requires R)
cat rms.table | rms2avglink.csh $rmscutoff > cluster
calc_score_for_cluster.pl cluster $scorefiles/$pdb.$topn.sc \
	  > clusterscores
sort -nr +1 clusterscores > clusterscores.bysize
pwd
head -20 clusterscores.bysize
#lp cluster.ps

# make links to top n clusters
i=1
while [ $i -le 9 ]; do
    f=$(head -$i clusterscores.bysize | tail -1 | gawk '{print $3}')
    echo linking cluster $i to $f
    ln -s $f cluster$i.pdb
    i=$[ i + 1 ]
done

# put the scores in a local file
cp $scorefiles/$pdb.$topn.sc ./
# make a scorefile with just the cluster centers
head -1 $pdb.$topn.sc > clusters.sc
i=1
while [ $i -le 9 ]; do
    f=$(head -$i clusterscores.bysize | tail -1 | gawk '{print $3}')
    grep $f $pdb.$topn.sc >> clusters.sc
    i=$[ i + 1 ]
done


cd $basedir




