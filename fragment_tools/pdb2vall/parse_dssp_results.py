#!/usr/bin/env python

from os import system, path
from sys import argv, exit, stderr
from optparse import OptionParser
#from get_fulllength_fasta import get_fulllength_fasta

''' In this script, we are going to parse out the lines of dssp_phi dssp_psi dssp_sa
    from the results of dssp2threestateSS.pl.'''


def parse_dssp_results( pdb_fn, arg_fasta_dict ):
    PDB2VALL_PATH = path.abspath(path.dirname(__file__)) + "/"
    if not PDB2VALL_PATH:
        stderr.write("ERROR: you should specify the path where your packages are first.\n")
        return 0

    print "def parse_dssp_results( pdb_fn, arg_fasta_dict ) has been called."
    missing_den_list = map( int, open( pdb_fn[:5] + ".fasta", "r").readline().split("missing_density_rsn:")[1:][0].strip().split())
    fasta_dict = arg_fasta_dict

    script_dssp2threestateSS = PDB2VALL_PATH + "structure_profile_scripts/dssp2threestateSS.pl"
    cmd = script_dssp2threestateSS + " " + pdb_fn + " > /dev/null 2> /dev/null"
    system( cmd )

    dssp_Dict  = {}

    dssp_results = open( pdb_fn + ".dssp", "r" )
    line = dssp_results.readline()
    count = 1

    while line:
        if "#  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA" in line:
            line = dssp_results.readline()
            while line:
                if not "!" in line:
                    if count in missing_den_list:
                        count += 1
                        continue
                    else:
                        newrsn = count
                        #rsn_in_dssp_results = int( line.strip().split()[1] )
                        rsd_in_dssp_results = line.strip().split()[3]

                        if newrsn > len( fasta_dict.keys()): break
                        if fasta_dict[ newrsn ] == rsd_in_dssp_results:
                            #if options.debug:
                            #    print count, fasta_dict[ newrsn ], rsd_in_dssp_results
                            dssp_sa  = line[35:38].strip()
                            dssp_phi = line[103:109].strip()
                            dssp_psi = line[109:115].strip()
                            #print newrsn, dssp_sa, dssp_phi, dssp_psi

                            dssp_Dict[ newrsn ] = [dssp_sa, dssp_phi, dssp_psi]

                            count += 1
                            line = dssp_results.readline()
                        else:
                            stderr.write("parse_dssp_results: ERROR: at pos %s pdbseqres( %s ) | dssp_rsd( %s )\n" % ( newrsn, fasta_dict[ newrsn ], rsd_in_dssp_results ))
                            exit()
                else:
                    line = dssp_results.readline()

        line = dssp_results.readline()

    return dssp_Dict


if __name__ == "__main__":
    parser = OptionParser()

    parser.add_option("-p", dest="pdb_fn", help="pdb file name eg 2oxgZ.pdb")
    parser.add_option("-d", action="store_true", dest="debug", help="print out some information")
    (options,args) = parser.parse_args()

    if not ( options.pdb_fn ):
        print "you miss pdb_fn, please look at info below"
        parser.print_help()
        exit()

    #parse_dssp_results( options.pdb_fn )
    #- change the arg back"
    print "this is a function being called by other script - should give the fasta_dict as one of the argument."

