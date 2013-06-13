#!/bin/bash

# antibody_repertoire_graft.bash repertoire_directory

target=model.pdb

if [ -z "$1" ]
then
    echo Usage: $(basename $0) repertoire_directory [force/quick]
    echo run antibody.py on a repertoire to create grafted structures
    echo force option will recaculate antibodies for which $target already exists
    echo quick option skips rosetta post-processing \(relax/idealize\)
    exit
fi

repdir=$1
if [ ! -d "$repdir" ]; then
	echo directory $repdir does not exist
	echo have you run antibody_repertoire.py?
	exit
fi
cd $repdir

if [ "$2" == "force" ]; then force=true; else force=false; fi
if [ "$2" == "quick" ]; then norelax="--quick=1"; fi

dirs=`ls -d *`
echo Processing Repertoire Antibodies: $dirs $norelax
success=0
failed=0
skipped=0

for d in $dirs; do
  if [ -d "$d" ]; then
  	echo "-------------------------------------- Antibody $d --------------------------------------"
    cd $d
 	if [ $force == true ] || [ ! -f grafting/$target ]
 	then 
	 	antibody.py --heavy-chain $d\H.fasta --light-chain $d\L.fasta --prefix grafting/ \
    	            $norelax 2>&1 | tee grafting.log
    	if [ $PIPESTATUS -eq 0 ]; then (( success++ )); else (( failed++ )); fi
	else
		echo $target exists, skipping antibody $d
		(( skipped++ ))
	fi
    cd ../
  fi
done

echo "-------------------------------------- Repertoire Complete --------------------------------------"
echo $(echo $dirs|wc|awk '{print $2}') Antibodies processed
echo $success Cases suceeded
echo $skipped Cases skipped
echo $failed Cases failed
