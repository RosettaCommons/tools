#!/usr/bin/env python

import sys
#from rosetta import rosetta,Pose,pose_from_pdb # no time improvement
from rosetta import *
from rosetta.protocols.antibody import *
rosetta.init(extra_options='-in:ignore_unrecognized_res')

def kink_begin(abinfo):
    H3_end = abinfo.get_CDR_loop(h3).stop()
    return H3_end - 2

def kink_end(abinfo):
    H3_end = abinfo.get_CDR_loop(h3).stop()
    return H3_end + 1

def kink_anion_res(abinfo):
    resi = abinfo.get_CDR_loop(h3).stop() - 1
    return resi

def kink_cation_res(abinfo):
    resi = abinfo.get_CDR_loop(h3).start() - 1
    return resi


def kink_constraints(pose, abinfo, outf = sys.stdout):
    kb=kink_begin(abinfo)
    kr0=kb
    kr1=kb+1
    kr2=kb+2
    kr3=kb+3
    kan=kink_anion_res(abinfo)
    kcat=kink_cation_res(abinfo)
    # CTER kink dihedral angle q
    outf.write("Dihedral CA %i CA %i CA %i CA %i SQUARE_WELL2 0.523 0.698 600\n" % (kr0,kr1,kr2,kr3) )
    # CTER kink q bond distance
    #outf.write("AtomPair CA %i CA %i FLAT_HARMONIC 7.125 0.5 0.625\n" % (kr0,kr3) )
    # KD Hbond
    outf.write("AtomPair N %i O %i FLAT_HARMONIC 2.0 2.0 2.0\n" % (kan,kcat) )
    return



def main(args):
    '''Output kink constraint files for a set of antibodies.  Usage: python kink_constraints.py [pdb_files]
    '''
    if len(args) < 2:
        print
        print main.__doc__
        return

    for filename in args[1:]:
        print filename

        try:
            pose = pose_from_pdb(filename)
            abinfo = AntibodyInfo(pose)
        except:
            print 'Problems reading %s, skipping.' % filename
            continue

        print pose
        print abinfo

        outfname = filename[:-4]+'.constr'
        outf = file(outfname, 'w')
        constraint = kink_constraints(pose,abinfo,outf)
        print "Constraints file %s:" % outfname
        outf.close()
        with open(outfname, 'r') as fin:
            print fin.read()

if __name__ == "__main__": main(sys.argv)

