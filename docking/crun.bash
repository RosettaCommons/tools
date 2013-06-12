if [ -z "$2" ]
then
    echo "Usage: crun.bash <pdb> <config>"
    echo "Sets up and launches a condor run"
    echo "<config> is the name of the condor config file"
    echo "without the '.config' extension"
    echo "The extension on the config file must be '.config'"
    exit
fi

pdb=$1

#Name of config file (not including extension)
#extension must be '.config'
config=$2

#Check for files before starting condor

if [ ! -f paths.txt ]
then
    echo paths.txt not found -- exiting!
    exit
fi

if [ ! -f $pdb.pdb ]
then
    echo $pdb.pdb not found -- exiting!
    exit
fi

if [ ! -f prepack/$pdb.ppk.pdb ]
then
    echo prepack/$pdb.ppk.pdb not found
    echo "You must prepack before running condor -- exiting"
    exit
fi

#Set up directories (replaces pdb_dir_maker.pl)
mkdir -p $pdb/outerr

exe=$pdb.$config.bash
cscript=$pdb.$config.con

#Get variables from the config file
source $config.config

command="rrun.sh $compiler $pdb $prefix $protocol_flags -nstruct $nstruct $sidechain_flags $search_flags $scorefilter $fab_flags $extra_flags"

#Set up an executable to run RosettaDock
echo \#\!/bin/bash > $exe

echo $command >> $exe
chmod +x $exe

#Make a relevant condor script
cp $rosetta_scripts/condor_scripts/base.con $cscript
echo "Executable = ./$exe" >> $cscript
echo "pdb = $pdb" >> $cscript
echo "Queue $Njobs" >> $cscript

condor_submit $cscript


