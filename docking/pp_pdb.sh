#!/bin/bash


if [ -z "$1" ]
then
    echo Usage: $(basename $0) pdb_code [push]
    echo post-process a pdb
    echo push flag transfers files back to hep
    exit
fi



PATH=$PATH:~/simcode/scripts

pdb=$1  

pp_compile_scorefiles.sh $pdb | tee $pdb.pplog
pp_extract_set.sh $pdb 10     | tee $pdb.pplog2
pp_extract_set.sh $pdb 200    | tee -a $pdb.pplog2
pp_extract_set.sh $pdb 1000    | tee -a $pdb.pplog2

if [ "$2" = "push" ]
then
  pp_push_set.sh $pdb
fi
