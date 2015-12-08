ROSETTA_PATH=/home/sid/simcode/rosetta_ensemble_final
SCRIPT_PATH=/home/sid/simcode/rosetta_scripts/docking

pdb=$1
chain_id=A

echo extracting out models from NMR structure
$SCRIPT_PATH/pdb_scripts/extract_model_from_pdb.pl $pdb

ls | grep "$pdb.m.*" > list1

cp $pdb.pdb $pdb.backup.pdb
cp $pdb.m.0.pdb $pdb.pdb

mkdir prepack
mkdir $pdb
cp paths.txt ./$pdb
cp $pdb.pdb ./prepack/
cp $pdb.*.pdb ./prepack/
mv list1 ./$pdb/

echo creating fasta file
$SCRIPT_PATH/pdb_scripts/pdb_fasta.pl $pdb.pdb > $pdb$chain_id.fasta

echo idealizing structures. . .
cd $pdb
$ROSETTA_PATH/rosetta.gcc ii $pdb $chain_id -idealize -nstruct 1 -l list1 

rm *.fasc
ls | grep "ii" > list2

#mv list2 ../prepack/
mv ii* ../prepack/

echo refining structures. . .
$ROSETTA_PATH/rosetta.gcc nn $pdb $chain_id -relax -farlx -skip_fragment_moves -nstruct 1 -l list2 

cp nn*.pdb ../
cd ..

rm -r prepack
rm -r $pdb


