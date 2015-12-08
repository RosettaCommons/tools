#!/bin/bash


if [ -z "$2" ]
then
    echo Usage: $(basename $0) \[-cf\] pdb_code top_set_size
    echo extract the top scoring structures
    echo   use -c to extract the *closest* structures to native
    echo   use -c to extract the top structures which are *far* from native
    exit
fi

#default to look at top scores
label=top
#masterext=top5000.sc
masterext=master

# load options
set -- $(getopt cf "$@")

while [ "$1" != "--" ]
do
    case "$1" in
	-c)  # get close structures
	    masterext=sc.byrms.top1000
	    label=close
	    shift;;
	-f)  # get close structures
	    masterext=sc.far1000
	    label=far
	    shift;;
    esac
done
shift 


basedir=$PWD
pdb=$1
ndecoys=$2
scorefiles=$basedir/$pdb/scorefiles

master=$pdb.$masterext
    


echo extracting $label $ndecoys structures of $pdb according to $master

nplus1=$[ $ndecoys + 1 ]
topn=$label$ndecoys
 
cd $scorefiles

head -$nplus1 $master > $pdb.$topn.sc
sort $pdb.$topn.sc -n +2 > $pdb.$topn.sc.byrms
head -$nplus1 $master |tail -$ndecoys | cut -c 1-20 > $pdb.$topn
adddirs.pl $pdb.$topn > $pdb.$topn.dirs

### make links
#cd ../$pdb
#lnset.sh $topn $(cat ../scorefiles/$pdb.$topn)

## make directory and unzip files
cd $basedir/$pdb
if [ -d $topn ]; then mv $topn $topn.$(date +%Y_%b%d_%H%M); fi
mkdir $topn
cd $topn
echo Unzipping files in $topn:
for f in $(cat $scorefiles/$pdb.$topn)
do
    prefix=$(echo $f | gawk '{print substr($1,0,2)}')
    # file could be in one of several places
    if [ -f $basedir/$pdb/$prefix$pdb.pdbs.tar.gz ] # zipped directory
      then tar zxvf $basedir/$pdb/$prefix$pdb.pdbs.tar.gz $f
    elif [ -f $basedir/$pdb/$prefix/$f ]            # in directory
      then cp $basedir/$pdb/$prefix/$f .
    elif [ -f $basedir/$pdb/$prefix.tgz ]           # zipped in directory
      then tar zxvf $basedir/$pdb/$prefix.tgz $prefix/$f;
	mv $prefix/$f $f; rmdir $prefix;
    else
      echo $f not found!!
    fi
done

cd $basedir
echo Extracted $(ls $pdb/$topn/|wc|cut -c1-8) of $ndecoys successfully



