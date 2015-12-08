#!/bin/bash


if [ -z "$1" ]
then
    echo Usage: $(basename $0) pdb_code
    echo push the \(small\) scorefiles and the top structures for
    echo a given pdb back to the main bast cluster
    exit
fi


basedir=$PWD
pdb=$1
scorefiles=$basedir/$pdb/scorefiles

#remote_basedir=/users/jeff/runs/$(basename $basedir)
remote_basedir=$basedir
remote_scorefiles=$remote_basedir/$pdb/scorefiles

ssh duke mv $remote_scorefiles $remote_scorefiles.$(date +%Y_%b%d_%H%M)
ssh duke mv $remote_basedir/$pdb $remote_basedir/$pdb.$(date +%Y_%b%d_%H%M)
ssh duke mkdir -p $remote_scorefiles
ssh duke mkdir -p $remote_basedir/$pdb

#  cd $scorefiles
#  keeper_exts="decoys.sc.sort master sc.byrms.top* top10* top200*"
#  keepers=$(for x in $keeper_exts; do echo $pdb.$x; done;)
#  scp -p $keepers duke.baker:$remote_scorefiles/
## keep everything now...in case of loss on remote system

cd $scorefiles
scp -p * duke:$remote_scorefiles/

cd $basedir
topdirs=$('ls' -d $pdb/top* | grep -v _ )
scp -rp $topdirs duke:$remote_basedir/$pdb/





