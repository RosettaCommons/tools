
if [ -z "$1" ]
then
    echo Usage: $(basename $0) pdb_code [extra_flags]
    echo Performs a standard prepacking of the monomer docking components
    echo including detection of disulfide bonds.  
    echo Starting pdb and paths.txt must be in the current directory, 
    echo and paths.txt should have ../ for starting structure path
    echo Extra_flags could be:
    echo \ \ -norepack1 or -norepack2
    echo \ \ -fab1      or -fab2
    exit
fi

pdb=$1
shift

if [ ! -f $pdb.pdb ]
then
    echo $pdb.pdb not found -- exiting!
    exit
fi

if [ ! -f paths.txt ]
then
    echo paths.txt not found -- exiting!
    exit
fi

#Insert -fab1 or -fab2 here if you are doing an antibody run.
extra_flags="$@"

mkdir -p prepack
cp $pdb.pdb prepack/

#If you make a new executable, e.g. rosetta.<compiler> change the variable
#$compiler to reflect this.
compiler=gcc  #default, new c++ version of rosetta

#Overall structure of docking protocol
pdb_flags="$pdb aa -s $pdb"

#Prepacking flags. used if docking is performed with rotamer trial minimization
protocol_flags="-prepack_rtmin -quiet"
#the commandline for prepacking structures for runs without rotamer trial minimization would be:
protocol_flags="-prepack_full -quiet"
#both options can be combined: this will invoke first a full repacking of all side chains in the complex, followed by a round of rotamer trial minimization 

#If you want to fix the sidechains of partner 1, use -norepack1
#If you want to fix the sidechains of partner 2, use -norepack2
sidechain_flags="-ex1 -ex2aro_only -unboundrot"

rrun.sh $compiler $pdb_flags $protocol_flags $sidechain_flags $extra_flags | tee $pdb.out

mv $pdb.out $pdb
cp $pdb/aa/$pdb.ppk.pdb prepack
mv $pdb $pdb.ppk

