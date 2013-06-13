#!/usr/bin/env bash

## prun.sh: same as rrun.sh without the ppk starting structure to rosetta, so we can run from the native for prepacking

# usage:
# rosettarun.sh compiler pdb_code

# makes a directory named for the pdb_code and runs in that directory.
# works with a Condor script or alone.

if [ $# -lt 3 ]
then
  echo ERROR--Not enough arguments: $@
  echo Usage: $0 compiler pdb_code prefix_or_number extra_flags
  exit 1
fi

InputLine=$@
COMPILER=$1
PDB=$2
prefix=$($rosetta_scripts/makeprefix.sh $3)
shift 3

#--- error checking -------------
if [ ! -f paths.txt ]; then echo No paths.txt ...exiting; exit;fi
if [ ! $prefix ]; then echo No prefix detected...exiting; exit;fi

#Disulfides (added by Mike Daily, 3/3/05)
#Calls makefixdisulf.py (found in $rosetta_scripts)
#Turns on disulfide flag if a fixdisulf is made, otherwise not
#Need to do before cd into $PDB
#====================================================================
makefixdisulf.py $PDB.pdb

#if [ -f $PDB.fixdisulf ]
#then
#    disulfflags="-fix_disulf ../$PDB.fixdisulf -use_disulf_logfile $PDB.disulflog -norepack_disulf"
#else
    disulfflags=""
#fi
#====================================================================

mkdir -p $PDB
cd $PDB

#--- prep pdb path---------------
# append prefix to pdb path
sed -e "s/^pdb path[[:print:]]*\//&$prefix\//" ../paths.txt > paths.txt
# prepare a directory for structures
pdbpath=$(grep 'pdb path' paths.txt | cut -c33-);mkdir -p $pdbpath


#--- Executable -----------------
exe="$rosetta_bin/rosetta.$COMPILER"
if [ ! -x $exe ]; then echo No executable $exe; exit; fi


#--- Arguments ------------------
# standard arguments plus optional arguments
args="$prefix $PDB _ -dock $@"
# add starting structure, unless it's in the optional arguments
if ! echo $args | grep -q "\-s "; then args="$args -s $PDB.ppk"; fi
# prepack mode should have at least chat-level output
# I have disabled this because -chat puts out a lot more than it used to
# Mike Daily, 3/3/05
#if echo $args | grep -q "\-prepack"; then args="$args -chat"; fi

args="$args $disulfflags"

echo --------------------------------
echo Rosetta Run in $PWD $(eval date)
echo $(eval uname -a)
echo $(basename $0) $InputLine
echo PDB: $PDB
echo prefix: $prefix
echo pdbpath: $pdbpath
echo exe: $exe
echo args: $args
echo --------------------------------


#--- Run ------------------------
fail="no"
if [ $CONDOR_VM ] && ! echo $args | grep -qE "\-verbose|\-inform|\-chat|\-yap"
then
    # Running inside condor, don't keep output unless -verbose flag is on
    /usr/bin/time nice $exe $args 2>/dev/null || fail=yes
else
    # Running interactively
    /usr/bin/time nice $exe $args 2>&1        || fail=yes
fi


echo -------------------
if [ "$fail" = "no" ]; then
echo Run finished \($PDB\)
else
echo Run failed!! \($PDB\)
fi
echo -------------------


#--- Clean -----------------------
# tar and move files from scratch directories
if echo $pdbpath | grep scratch
then
    savedir=$PWD
    echo pdb files in $pdbpath are being zipped and moved to $savedir
    cd $pdbpath
    tarfile=$prefix$PDB.pdbs.tar.gz
    tar -czvf $tarfile *$prefix*pdb --remove-files
    mv $tarfile $savedir
fi

cd ../

exit
