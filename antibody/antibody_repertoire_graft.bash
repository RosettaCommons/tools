#!/bin/bash

# antibody_repertoire_graft.bash repertoire_directory

target=model.pdb

if [ -z "$1" ]
then
    echo Usage: $(basename $0) repertoire_directory [force]
    echo run antibody.py on a repertoire to create grafted structures
    echo force option will recaculate antibodies for which $target already exists
    exit
fi

repdir=$1
cd $repdir
dirs=`ls -d *`
echo Processing Repertoire Antibodies: $dirs
success=0
failed=0
skipped=0

if [ "$2" == "force" ]; then force=true; else force=false; fi

for d in $dirs; do
  if [ -d "$d" ]; then
  	echo "-------------------------------------- Antibody $d --------------------------------------"
    cd $d
 	if [ $force == true ] || [ ! -f grafting/$target ]
 	then 
	 	antibody.py --heavy-chain $d\H.fasta --light-chain $d\L.fasta --prefix grafting/ \
    	            --relax=1 --idealize=0 \
        	        2>&1 | tee grafting.out
        	        # --rosetta-platform linuxgccrelease \
    	if [ $PIPESTATUS -eq 0 ]; then (( success++ )); else (( failed++ )); fi
#    	ln -sF grafting/model.pdb $d.pdb
#    	ln -sF grafting/grafted.pdb $d.grafted.pdb
#    	ln -sF grafting/grafted.relaxed.pdb $d.grafted.relaxed.pdb
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
