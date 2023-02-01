#!/usr/bin/env python
# quick and dirty PDB reader

def read_pdb( filename ):

    coords = {}
    pdb_lines = {}
    sequence = {}

    old_resnum = 0
    old_chain  = ''
    old_segid = '    '
    chains = []
    residues = []
    segids = []
    for line in open( filename ):

        if (len(line)>54 and  (line[0:4] == 'ATOM' or line[0:4] == 'HETA' ) ):

            resnum = int( line[22:26] )
            chain = line[21]
            segid = line[72:76]
            atom_name = line[12:16]
            position = [float(line[30:38]),float(line[38:46]),float(line[46:54])]

            if not ( chain in coords.keys() ):
                coords[chain] = {}
                pdb_lines[chain] = {}
                sequence[ chain ] = {}
            if not (segid in coords[chain].keys()):
                coords[chain][segid] = {}
                pdb_lines[chain][segid] = {}
                sequence[ chain ][segid] = {}

            sequence[chain][segid][resnum] = line[17:20]

            if not ( resnum in coords[chain][segid].keys() ):
                coords[chain][segid][resnum] = {}
                pdb_lines[chain][segid][resnum] = {}

            coords[chain][segid][resnum][atom_name] = position
            pdb_lines[chain][segid][resnum][atom_name] = line[:-1]

            if ( len(residues) == 0 or \
                 resnum != old_resnum or chain != old_chain or segid != old_segid ):
                chains.append( chain )
                residues.append( resnum )
                segids.append( segid )
            old_resnum = resnum
            old_chain  = chain
            old_segid  = segid

    return ( coords, pdb_lines, sequence, chains, residues, segids )
