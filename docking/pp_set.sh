#!/bin/bash

# post-process all pdbs in the current directory

pdbs=$(ls -d ???? | grep "^[1-9|c]")
echo PDB set for postprocessing: $pdbs
if [ "$@" = "push" ]; then echo Pushing data to bast; fi

#getAllscores.pl > Allscores

for pdb in $pdbs
do
    echo Postprocessing $pdb
    pp_pdb.sh $pdb $@
done
