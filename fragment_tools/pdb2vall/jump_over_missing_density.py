#!/usr/bin/env python

from sys import argv, exit, stderr
from os.path import exists

def jump_over_missing_density( fasta_fn, pdb_fn, request ):
    missing_density_rsn_list = open( fasta_fn, "r").readline().split("missing_density_rsn:")[1:][0].strip().split()
    pdb_xyz_file     = []
    pdb_torsion_file = []

    ## RENUMBER IDEALIZED PDB BACK TO PDB_SEQRES
    pdb_fn = pdb_fn[:5]
    if exists( pdb_fn + '_0001_0001.pdb' ):

        file = open( pdb_fn + '_0001_0001.pdb', "r")
        line = file.readline()

        newnum_xyz     = 1
        newnum_torsion = 1

        while line:
            if line[12:16] == " CA ":
                if str(newnum_xyz) in missing_density_rsn_list:
                    newnum_xyz = newnum_xyz + 1
                    continue

                xcord = line[30:38]
                ycord = line[38:46]
                zcord = line[46:54]
                xyz_line_edit = "%s%4d %s" %( line[:22], newnum_xyz, line[27:] )
                pdb_xyz_file.append( xyz_line_edit )
                newnum_xyz = newnum_xyz + 1

            if line.startswith("REMARK") and "torsion" not in line:
                if str(newnum_torsion) in missing_density_rsn_list:
                    newnum_torsion = newnum_torsion + 1
                    continue

                idealized_rsd = line.split()[4]
                if idealized_rsd != "X":
                    secstr        = line.split()[5]
                    phi           = line.split()[6]
                    psi           = line.split()[7]
                    omega         = line.split()[8]

                    torsion_line_edit = "%s %s %s %s %s %s\n" %( newnum_torsion, idealized_rsd, secstr, phi, psi, omega )
                    pdb_torsion_file.append( torsion_line_edit )
                    newnum_torsion = newnum_torsion + 1

            line = file.readline()
        file.close()
    else:
        stderr.write("ERROR: there's no idealized and relaxed %s_0001_0001.pdb here!\n"  % pdb_fn )
        return 0
        exit()


    if request == "xyz":
        return pdb_xyz_file
    if request == "torsion":
        return pdb_torsion_file

if __name__ == '__main__':
    if len( argv ) < 3:
        print
        print "USAGE: %s fasta_fn pdb_fn xyz/torsion" % argv[0]
        print
        print "-"*75
        exit()

    print jump_over_missing_density( argv[1], argv[2], argv[3] )


