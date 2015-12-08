#!/bin/bash

if [ -z "$1" ]
then
    echo Usage: $0 pdb_code
    echo collect the scorefiles and sort them various ways
    echo execute from the root work directory \(above pdb_code\/\)
    exit
fi

basedir=$PWD
pdb=$1
scorefiles=$basedir/$pdb/scorefiles

cd $pdb || { echo directory $pdb not found...exiting; exit; }

## make scorefiles directory
if [ -d $scorefiles ]; then mv $scorefiles $scorefiles.$(date +%Y_%b%d_%H%M); fi
mkdir -p $scorefiles

## be sure temp directory exists for sorting
#if [ $TMPDIR -a ! -x $TMPDIR ]; then mkdir -p $TMPDIR; fi

## merge all score files, place in scorefiles directory
echo Merging $(ls *sc|wc|cut -c1-9) scorefiles
head -1 $(ls *sc |head -1) > $scorefiles/$pdb.sc
nfields=$(wc $scorefiles/$pdb.sc | gawk '{print $2;}')
cat *sc | \
    column_filter.pl $nfields | \
    tr -cd "[:print:]|[:space:]" | grep -v ZZ | \
    grep -v filename | grep -v '\*\*\*' | grep -v NaN >> $scorefiles/$pdb.sc
wc $scorefiles/$pdb.sc

## extract native & start from one run
echo -n Extracting native and input...
head -7 $(ls *sc | head -1) > $scorefiles/$pdb.native.sc
echo done

cd $scorefiles

## sort by rms
if grep -q native $scorefiles/$pdb.native.sc
then
  echo -n Sorting by rms...
  head -1 $pdb.sc > $pdb.sc.byrms
  filter_on.pl rms lt 12 $pdb.sc | sort -n +2 | grep -v filename >> $pdb.sc.byrms
  grep -v nat $pdb.sc.byrms | grep -v inp | head -1000 > $pdb.sc.byrms.top1000
  echo done
  echo Closest:
  head $pdb.sc.byrms.top1000 | cut -c 1-39,317-
else
  echo no native present, skipping rms sort
fi

## get decoys
echo -n Extracting decoys...
head -1 $pdb.sc > $pdb.decoys.sc
grep decoy $pdb.sc | grep -v NaN >> $pdb.decoys.sc
echo done
wc $pdb.decoys.sc

## sort
head -1 $pdb.sc > $pdb.decoys.sc.sort
sort -n +1 $pdb.decoys.sc >> $pdb.decoys.sc.sort

## get a top subset
ndecoys=5000
nplus1=$[ $ndecoys + 1 ]
topn=top$ndecoys
head -$nplus1 $pdb.decoys.sc.sort > $pdb.$topn.sc

## debump
filter_bumps.pl $pdb.$topn.sc > $pdb.$topn.nobumps.sc
wc $pdb.$topn.nobumps.sc

## and use this subset in future scripts
## 7/02: maybe bump filter is unnecessary these days...
rm -f $pdb.master
#ln -s $pdb.$topn.nobumps.sc $pdb.master # bump filter is lousy 2/14/02
ln -s $pdb.$topn.sc $pdb.master
ls -l $pdb.master
echo Best scores:
head $pdb.master | cut -c 1-39,317-

## get decoys far from native
if grep -q native $scorefiles/$pdb.native.sc
then
    filter_on.pl rms gt 12 $pdb.master | head -1000 > $pdb.sc.far1000
    echo done
fi
