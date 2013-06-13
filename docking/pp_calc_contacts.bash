#!/bin/bash

# assuming the top200 directories have Nmodels in pdbs labeled clusterX.pdb, 
# calculate the contacts and get the rmsds into the clustersummary file

if [ -z "$1" ]
then
    echo Usage: $(basename $0) pdb_code
    echo calculate the correct contacts
    echo currently hard-wired for the Chen benchmark set
    exit
fi

base_dir=$PWD
pdb=$1
topn=top200

# Chen benchmark file locations
native_dir=~/pdb/chen/co_prepped
map_dir=~/pdb/chen/maps
native_contacts_dir=~/pdb/chen/native_contacts

Nmodels=25

echo Calculating the correct contacts on the top $Nmodels of $pdb in $base_dir \($topn\)

#pdbs=$('ls' -d ????)
#echo PDBS: $pdbs
#
#for pdb in $pdbs
#  do
##### pdb given
  
  map=$map_dir/$pdb.map
  native=$native_dir/$pdb\_c.pdb
  
  cd $base_dir/$pdb/$topn
# loop over models -- find interfaces
  n=1
  while [ $n -le $Nmodels ]
    do
    model=cluster$n.pdb
    echo $model
    find_contacts.pl $model  > $model.contacts
    map_residues.pl $map $model.contacts > $model.contacts.mapped
    n=$[ $n + 1 ]
  done
  
# loop over models -- count
  rm -f clustercontacts
  n=1
  while [ $n -le $Nmodels ]
    do
      model=cluster$n.pdb
      count_contact_matches.pl $model.contacts.mapped $native.contacts > tmp1
      count_contact_matches.pl $native.contacts.mapped $model.contacts > tmp2
      
      echo "$model: $(cat tmp1) model, $(cat tmp2) native" | tee -a clustercontacts
    n=$[ $n + 1 ]
      
  done
  
  head -$Nmodels clusterscores.bysize > tmp1
  echo $pdb: > clustersummary
  paste tmp1 clustercontacts >> clustersummary
  echo --------------------------------------- >> clustersummary
  cat clustersummary

  rm tmp* -f
  cd ../../

#done

cd $base_dir
cat ????/$topn/clustersummary > summary.txt
