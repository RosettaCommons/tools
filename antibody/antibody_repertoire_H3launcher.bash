#!/bin/bash

# antibody_repertoire_H3_launcher.bash repertoire_directory

if [ -z "$1" -o "help" == "$1" -o "-h" == "$1" -o "-help" == "$1" -o "--help" == "$1" ]
then
    cat << EOUSAGE
NAME

$(basename $0) - perform H3 modelling for a repertoire of antibodies

SYNOPSIS

$(basename $0) repertoire_directory [force]

DESCRIPTION

The ab initio modelling of the H3 loop is particularly compute
intensive. This is implemented by the tool 'antibody_H3' of
Rosetta and may better run on a high-performance compute cluster
than on a local machine when performing on a large set of sequences.

The script is prepared for the developers' compute environment
but also anybody else with a ready installation of the queuing system
'slurm' experiences an immediate benefit.

SEE ALSO

 * slurm
 * antibody_H3
 * abH3.(qsub|sbatch)
 * abH3.flags

EOUSAGE
    exit -1
fi

repdir=$1
if [ "$2" == "force" ]; then force=true; else force=false; fi
maxjobs=200

cd $repdir
if [[ -f PDBlist ]]
then
	dirs=`cat PDBlist`
else
	dirs=`ls -dp * | grep '/' | sed 's/\///'`
fi

echo Processing Repertoire Antibodies: $dirs
complete=0
queued=0
launched=0

launchscript=abH3.qsub
if [[ `hostname` = *stampede* ]]
then
	launchscript=abH3.sbatch
fi

if [ ! -f abH3.flags ] || [ ! -f $launchscript ]
then
	echo Missing abH3.flags or $launchscript, exiting...
	exit 1
fi

lastdecoynum=`grep -E '^-nstruct' abH3.flags | awk {'print $2'}`
lastdecoyfile=pdbs/model_$lastdecoynum.pdb
echo Seeking $lastdecoynum decoys for each antibody
echo Allowing a maximum of $maxjobs jobs in the queue

export SQUEUE_FORMAT="%.7i %.9P %.21j %.8u %.2t %.10M %.6D %R"


for d in $dirs; do
	if [ -d "$d" ]; then
  		echo -n "Antibody $d..."
    	cd $d
    	if [ -f $lastdecoyfile ] || [ -f $lastdecoyfile.gz ]
    	then
			echo $lastdecoyfile exists, complete
			(( complete++ ))
			cd ../
			continue
		fi
		jobname=$repdir-$d
    	squeue -u $USER | grep -q $jobname
    	jobinqueue=$?
    	if [ $jobinqueue = 0 ]
    	then
    		echo in queue, skipping
    		(( queued++ ))
    	else
 	  		echo
 	  		if [ ! -d outerr ]; then mkdir outerr; fi
 	  		if [ ! -d pdbs ]; then mkdir pdbs; fi
 	   		sed 's/ABNAME0000/'$jobname'/' ../$launchscript > $launchscript
 	    	sbatch $launchscript
 	    	echo
 	    	(( queued++ ))
 	    	(( launched++ ))
 	    fi
		cd ../
    	if [ $queued = $maxjobs ]
    	then
    		echo maximum jobs submitted \($maxjobs\)
	    	echo Exiting
    		break
    	fi
  	fi
done

total=$(echo $dirs|wc|awk '{print $2}')
(( abovemax = $total - $complete - $queued ))

echo "-------------------------------------- Repertoire Complete --------------------------------------"
echo $total Antibodies processed
echo $complete Jobs complete
echo $launched Jobs newly launched
echo $queued Jobs enqueued
echo $abovemax Jobs remain unlaunched
echo "-------------------------------------------------------------------------------------------------"
squeue -u $USER


##Restarting with power!
#cd ../3GI9
#sed 's/ABNAME0000/'$(basename $pwd)'/' ../abH3.power.sbatch > abH3.power.sbatch
#sbatch abH3.power.sbatch

