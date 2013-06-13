#!/bin/sh

pdb=$1
voidoo=~jeff/bin/voidoo

# this line removes the .pdb from the filename, if there is one
pdb=$(echo $pdb | sed 's/\.pdb//')
  
sed "s/PDB/$pdb/g" ~jeff/etc/voidoo.script | $voidoo
grep calcns $pdb.log | sort -nr +6 > $pdb.voids
