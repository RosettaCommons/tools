#!/bin/bash

# antibody_repertoire_graft.bash repertoire_directory

target=model.pdb

if [ -z "$1" -o "-h" == "$1" -o "help" == "$1" -o "-help" == "$1" -o "--help" == "$1" ]
then
    cat <<EOUSAGE

USAGE: $(basename $0) repertoire_directory [force/quick]

DESCRIPTION

$(basename $0) runs antibody.py on every single of an arbitrarily
large collection of antibody seqeuences, referred to as 'repertoire',
to create grafted structures. The repertoire is expected to be represented
as a directory in which every entry is another directory with the name
of the antibody. This in turn needs to contain the two FASTA files
'dirnameL.fasta' and 'dirnameH.fasta' that describe the light and heavy
chain of the antibody. The script 'antibody_repertoire.py' helps preparing
the setup.

The first argument specifies the directory representing the repertoire.
The second/third argument may be 'force' and/or 'quick':

 * force will recaculate antibodies for which the file $target already exists
 * quick skips rosetta post-processing (relax/idealize)

OUTPUT

The tool antibody.py is run within the directory harboring the information
for a particular antibody. The output goes into the subdirectory 'grafting'
and receives the file name '$target'.

SEE ALSO

 * antibody.py
 * antibody_repertoire.py

EOUSAGE
    exit -1
fi

repdir=$1
if [ ! -d "$repdir" ]; then
	echo directory $repdir does not exist
	echo have you run antibody_repertoire.py?
	exit
fi
cd $repdir

if [ "$2" == "force" -o "$3" == "force"]; then force=true; else force=false; fi
if [ "$2" == "quick" -o "$3" == "quick"]; then norelax="--quick=1"; fi

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
