#!/bin/bash


if [ -z "$1" ]
then
    echo Usage: $(basename $0) \[-extrafilter score gt\/lt val\] \[-bktot\] pdb_code_or_ALL \[rms=5\]
    echo examine the dwindling of close \(\<rms\) structures through the filters
    echo -extrafilter flag indicates whether to filter for an extra score \>\/\< x
    echo -bktot flag indicates a sort by bk_tot score
    exit
fi

extrafilter=-999
extramsg=
extravar=contact
if [ "$1" = "-extrafilter" ] # sort on bktot
then   
    extravar=$2
    if [ "$3" = "lt" ]
    then 
	extrasym="<"
	extraltgt=$3
    else
        extrasym=">"
	extraltgt="gt"
    fi
    extrafilter=$4
    extramsg=" with $extravar$extrasym$extrafilter"
    shift 4
fi

# default: sort on score
sorton=score
sortext=byscore
sortcol=1
if [ "$1" = "-bktot" ] # sort on bktot
then   
    sorton=bk_tot
    sortext=bybktot
    sortcol=23
    shift
fi

if [ -z "$2" ]
then
    rms=5
else
    rms=$2
fi
    
basedir=$PWD

if [ "$1" = "ALL" ]
then
    pdbs=$('ls' -d ????)
else
    pdbs=$1
fi

echo Sorting on $sorton$extramsg\; Number of decoys within $rms Angstroms in $pdbs:
printf "%10s %10s %10s %10s %10s %10s %10s\n" "pdb" "basis" "all" "40000" "4000" "200" "10"

for pdb in $pdbs
do

    scorefiles=$basedir/$pdb/scorefiles
    if [ ! -x $scorefiles ]; then scorefiles=$basedir/scorefiles; fi
    if [ ! -x $scorefiles ]
    then
	printf "%10s scorefiles not found\n" $pdb
    else

	cd $scorefiles

	printf "%10s" $pdb

	# make sorted file, if necessary
	masterf=$pdb.sc.$sortext
	if [ ! -f $masterf ]
	then
	    #echo Sorting by score to $masterf
	    head -1 $pdb.sc > $masterf
	    grep -vE 'filename|nat|inp' $pdb.sc | sort +$sortcol -n >> $masterf
	fi

	basis=$(filter_on.pl $extravar $extraltgt $extrafilter $masterf 2>/dev/null | \
		grep -v 'filename' | wc | gawk {'print $1;'})
	printf " %10d" $basis

	# count the hits
	for nt in 40000001 40001 4001 201 11
	do 
	    n=$(head -$nt $masterf | \
		filter_on.pl rms lt $rms 2>/dev/null | \
		filter_on.pl $extravar $extraltgt $extrafilter 2>/dev/null | \
		grep -v 'filename' | wc | gawk {'print $1;'})
	    printf " %10d" $n
	done
	
	printf "\n"
	cd $basedir

    fi

done




#  cd runs/non-10-29
#  pull_everywhere "runs/non-10-29/dwindle*" .
#  cd ../HA-10-25
#  pull_everywhere "runs/non-10-29/dwindle*" .
#  cp ../non-10-29-HSC/dwindle* .

#  dwindle=dwindlebk
#  cat $dwindle.* > tmp; head -2 tmp > $dwindle-ALL; grep -vE 'Number|found|pdb' tmp |sort >> $dwindle-ALL; echo >> $dwindle-ALL; grep -vE 'Number|found|pdb' ../HA-10-25/$dwindle.* |sort >> $dwindle-ALL; cat $dwindle-ALL



