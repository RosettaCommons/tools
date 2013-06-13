ROSETTA_PATH=/home/sid/simcode/rosetta_ensemble_final
SCRIPT_PATH=/home/sid/simcode/rosetta_scripts/docking

pdb=$1
chain_id=$2
n_struct=$3

echo building $n_struct conformers from $pdb.pdb $2

cp $pdb.pdb $pdb.backup.pdb
$SCRIPT_PATH/pdb_scripts/pdb_fasta.pl $pdb.pdb > $pdb$chain_id.fasta

mkdir $pdb
mkdir prepack
cp $pdb.pdb ./prepack/
cp paths.txt ./$pdb/

cd $pdb
$ROSETTA_PATH/rosetta.gcc ii $pdb $chain_id -idealize -nstruct 1 -s $pdb.pdb
inum=_0001

cp ii$pdb$inum.pdb ../prepack/$pdb.pdb

$ROSETTA_PATH/rosetta.gcc rr $pdb $chain_id -relax -skip_fragment_moves -nstruct $n_struct -s $pdb.pdb

cp rr*.pdb ../
cd ..
rm -r $pdb
rm -r prepack 
rm $pdb.backup.pdb
rm $pdb$chain_id.fasta
