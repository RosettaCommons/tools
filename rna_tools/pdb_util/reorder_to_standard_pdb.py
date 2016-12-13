#!/usr/bin/env python
from read_pdb import read_pdb
from sys import argv
pdb = argv[1]

( coords, pdb_lines, sequence, chains, residues ) = read_pdb( pdb )


atoms = {
        "  A":[" P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," O2'"," C1'"," N9 "," C8 "," N7 "," C5 "," C6 "," N6 "," N1 "," C2 "," N3 "," C4 "],
        "  G":[" P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," O2'"," C1'"," N9 "," C8 "," N7 "," C5 "," C6 "," O6 "," N1 "," C2 "," N2 "," N3 "," C4 "],
        "  U":[" P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," O2'"," C1'"," N1 "," C2 "," O2 "," N3 "," C4 "," O4 "," C5 "," C6 " ],
        "  C":[" P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," O2'"," C1'"," N1 "," C2 "," O2 "," N3 "," C4 "," N4 "," C5 "," C6 "],
}


for (chain,residue) in zip(chains, residues):
        atom_keys = pdb_lines[ chain ][ residue ].keys()
        res = pdb_lines[ chain ][ residue ][ atom_keys[0] ][17:20]
        assert( res in atoms.keys() )
        for atom in atoms[ res ]:
                if ( not atom in atom_keys ):
                        print atom, "missing from", chain, residue
                        exit( 0 )
                print pdb_lines[ chain ][ residue ][ atom ]
