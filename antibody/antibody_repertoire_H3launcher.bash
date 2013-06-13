#!/bin/bash

# antibody_repertoire_H3_launcher.bash repertoire_directory

if [ -z "$1" ]
then
    echo Usage: $(basename $0) repertoire_directory [force]
    echo launch jobs for H3 modeling
    exit
fi

repdir=$1
if [ "$2" == "force" ]; then force=true; else force=false; fi
maxjobs=20

cd $repdir
dirs=`ls -dp * | grep '/' | sed 's/\///'`
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

