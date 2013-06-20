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

def trp_residue(abinfo):
    Wi = abinfo.get_CDR_loop(h3).stop() + 1
    return Wi

def kink_constraints(pose, abinfo, cst = sys.stdout, cst_fa = sys.stdout):
    # kink residues
    kb=kink_begin(abinfo)
    kr0=kb
    kr1=kb+1
    kr2=kb+2
    kr3=kb+3
    # CTER kink dihedral angle q
    cst.write(   "Dihedral CA %i CA %i CA %i CA %i SQUARE_WELL2 30 40 600 DEGREES\n" % (kr0,kr1,kr2,kr3) )
    cst_fa.write("Dihedral CA %i CA %i CA %i CA %i SQUARE_WELL2 30 40 600 DEGREES\n" % (kr0,kr1,kr2,kr3) )
    # CTER kink q bond distance
    #cst.write("AtomPair CA %i CA %i FLAT_HARMONIC 7.125 0.5 0.625\n" % (kr0,kr3) )

    # KD Hbond
    kan=kink_anion_res(abinfo)
    kcat=kink_cation_res(abinfo)
    cst.write(   "AtomPair N %i O %i FLAT_HARMONIC 2.0 2.0 2.0\n" % (kan,kcat) )
    cst_fa.write("AtomPair N %i O %i FLAT_HARMONIC 2.0 2.0 2.0\n" % (kan,kcat) )

    # KD sc Hbond
    if (pose.residue(kan).name3() == "ASP" or pose.residue(kan).name3() == "GLU") and (pose.residue(kcat).name3() == "ARG" or pose.residue(kcat).name3() == "LYS"):
        if pose.residue(kan).name3() == "ASP":
            kan_at = "OD1"
        else:  # GLU
            kan_at = "OE1"
        if pose.residue(kcat).name3() == "ARG":
            kcat_at = "NH1"
        else:  # LYS
            kcat_at = "NZ"
        cst_fa.write("AtomPair %s %i %s %i FLAT_HARMONIC 2.0 2.0 2.0\n" % (kan_at,kan,kcat_at,kcat) )

    # Trp Hbond
    Wi = trp_residue(abinfo)
    W  = pose.residue(Wi)
    if W.name3() == "TRP":
        cst_fa.write("AtomPair NE1 %i O %i FLAT_HARMONIC 2.0 2.0 2.0\n" % (Wi,kr1) )

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

        cstname = filename[:-4]+'.cst'
        cst = file(cstname, 'w')
        cst_faname = filename[:-4]+'.cst_fa'
        cst_fa = file(cst_faname, 'w')
        constraint = kink_constraints(pose,abinfo,cst,cst_fa)
        print "Constraints files %s, %s:" % (cstname, cst_faname)
        cst.close()
        cst_fa.close()
        with open(cst_faname, 'r') as fin:
            print fin.read()

if __name__ == "__main__": main(sys.argv)

