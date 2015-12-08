if [ -z "$1" ]
then
    echo Usage: $(basename $0) pdb_code [extra_flags]
    echo Performs a standard full atom run of the monomer docking components
    echo Starting pdb and paths.txt must be
    echo in the current directory.
    echo starting structure ppk.pdb in prepack
    echo Extra_flags could be:
    echo \ \ -norepack1 or -norepack2
    echo \ \ -fab1      or -fab2
    echo search_flags could be
    echo -randomize1, -randomize2
    echo or -dock_pert 3 8 8, -spin, etc.
    exit
fi

pdb=$1
shift

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

#Insert -fab1 or -fab2 here if you are doing an antibody run.
extra_flags="$@"

#If you make a new executable, e.g. rosetta.<compiler> change the variable
#$compiler to reflect this.
compiler=gcc  #default, new c++ version of rosetta

#Overall structure of docking protocol
pdb_flags="$pdb aa"

#Remove -dock_rtmin if you don't want to do rotamer minimization
#Change nstruct if you want to do more than one structure.
protocol_flags="-dock_mcm -dock_rtmin -quiet -nstruct 2000 -fake_native"

#How to do the docking search.  You can change the numbers on -dock_pert
#or use -randomize1 or -randomize2 instead of -dock_pert.
search_flags="-dock_pert 5 15 15 -spin"

#If you want to fix the sidechains of partner 1, use -norepack1
#If you want to fix the sidechains of partner 2, use -norepack2
sidechain_flags="-ex1 -ex2aro_only -unboundrot"

rrun.sh $compiler $pdb_flags $protocol_flags $search_flags $sidechain_flags $disulfflags $extra_flags | tee $pdb.out

mv $pdb.out $pdb
