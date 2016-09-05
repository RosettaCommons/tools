#!/usr/bin/env python

from sys import argv, exit, stderr
from os.path import exists
from numbering_back_to_pdbseqres import numbering_back_to_pdbseqres
from optparse import OptionParser


def renumber_structure_profile_checkpoint_back_to_seqres( fasta_fn, structure_profile_checkpoint_fn, arg_fasta_dict ):
    fasta_dict = arg_fasta_dict

    missing_den_list = map( int, open( fasta_fn, "r").readline().split("missing_density_rsn:")[1:][0].strip().split())
    checkpoint_count = 1
    checkpoint_dict  = {}
    file = open( structure_profile_checkpoint_fn, "r")
    line = file.readline()

    while line:
        if len(line) > 10:
            if checkpoint_count in missing_den_list:
                checkpoint_count += 1
                continue
            else:
                newrsn = checkpoint_count
                line_edit = line.strip().split()[1:]
                #for the reason that some jump over numbering
                if newrsn > len( fasta_dict.keys()):
                    break
                if fasta_dict[ newrsn ] == line.strip().split()[0]:
                    '''if options.debug:
                        print newrsn, line_edit
                        print "rsd at %s in seqres_fasta -->" % newrsn, fasta_dict[ newrsn ], line.strip().split()[0], "<-- rsd at %s in checkpoint_file" % newrsn
                        print "fasta_from_seqres %s ;checkpoint %s\n" % ( fasta_dict[ newrsn ], line.strip().split()[0] )'''
                    checkpoint_dict[ newrsn ] = " ".join( "%4.3f" % float(x) for x in line_edit)
                    checkpoint_count += 1
                else:
                    stderr.write( "ERROR: checkpoint rsd doesn't match the fasta from seqres, it stops at rsn %s\n" % newrsn )
                    stderr.write( "in which fasta_from_seqres thinks it should be %s while checkpoint thinks it should be %s\n" % ( fasta_dict[ newrsn ], line.strip().split()[0] ))
                    return 0
                    exit()

        line = file.readline()
    file.close()

    return checkpoint_dict


if __name__ == "__main__":
    parser = OptionParser()

    parser.add_option("-f", dest="fasta_fn", help="fasta file name")
    parser.add_option("-c", dest="structure_profile_checkpoint_fn", help="checkpoint file generated from structure profile")
    parser.add_option("-d", action="store_true", dest="debug", help="print out some information")
    (options,args) = parser.parse_args()

    if not ( options.fasta_fn or options.structure_profile_checkpoint_fn ):
        print "you miss some files, please look at info below"
        parser.print_help()
        exit()

    print "remember to change stuff back"
    #print renumber_structure_profile_checkpoint_back_to_seqres( options.fasta_fn, options.structure_profile_checkpoint_fn )
