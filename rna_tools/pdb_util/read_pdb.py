#!/usr/bin/env python
# quick and dirty PDB reader

def read_pdb( filename ):

    coords = {}
    pdb_lines = {}
    sequence = {}

    old_resnum = 0
    old_chain  = ''
    chains = []
    residues = []
    for line in open( filename ):

        if (len(line)>54 and  (line[0:4] == 'ATOM' or line[0:4] == 'HETA' ) ):

            resnum = int( line[22:26] )
            chain = line[21]
            atom_name = line[12:16]
            position = [float(line[30:38]),float(line[38:46]),float(line[46:54])]

            if not ( chain in coords.keys() ):
                coords[chain] = {}
                pdb_lines[chain] = {}
                sequence[ chain ] = {}

            sequence[chain][resnum] = line[17:20]

            if not ( resnum in coords[chain].keys() ):
                coords[chain][resnum] = {}
                pdb_lines[chain][resnum] = {}

            coords[chain][resnum][atom_name] = position
            pdb_lines[chain][resnum][atom_name] = line[:-1]

            if ( len(residues) == 0 or \
                 resnum != old_resnum or chain != old_chain ):
                chains.append( chain )
                residues.append( resnum )
            old_resnum = resnum
            old_chain  = chain

    return ( coords, pdb_lines, sequence, chains, residues )
