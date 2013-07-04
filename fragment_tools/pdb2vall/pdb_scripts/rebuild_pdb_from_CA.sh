#!/bin/bash
Bin=$(dirname $(readlink -f $0))
cwd=$(pwd)
pdb=$(basename $1 .pdb )

cp $1 $Bin/maxsprout/$pdb.alpha
cd $Bin/maxsprout/test/
runmaxsprout.pl ../$pdb.alpha $pdb.rebuilt
cd $cwd
mv  $Bin/maxsprout/test/$pdb.rebuilt ./$pdb.pdb.rebuilt
rm $Bin/maxsprout/$pdb.alpha
rm $Bin/maxsprout/test/$pdb*.log
