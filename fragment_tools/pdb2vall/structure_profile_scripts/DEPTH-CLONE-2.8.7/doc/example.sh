# compute residue depth for 1MVT.pdb, 10 iterations, cavity-size = 3 solvents, keep minimally solvated structures
./DEPTH -i 1MVT.pdb -o 1MVT -n 10 -survive 3 -keep 1MVT-sol

# compute change of residue depth upon chains forming complexes
./complex-depth -i 1LXY.pdb -o 1LXY.depth_dif

# Calculate residue solvent accessibility, 92 points on sphere, 1.4 Angstrom probe size
./ASA 1MVT.pdb 92 1.4 1MVT.asa

# predict-binding-site of 1MVT.pdb, binding cavity defined with 3 solvents.
./predict-binding-site -d 1MVT-residue.depth -a 1MVT.asa -p par/:3 -o 1MVT.cavity3.prediction -c 1MVT.pdb -y 1MVT.cavity3.pdb
