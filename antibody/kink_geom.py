#!/usr/bin/env python

import math
import sys
from rosetta import *
from rosetta.protocols.antibody import *
rosetta.init(extra_options='-in:ignore_unrecognized_res')

def kink_begin(abinfo):
    H3_end = abinfo.get_CDR_loop(h3).stop()
    return H3_end - 2

def kink_end(abinfo):
    H3_end = abinfo.get_CDR_loop(h3).stop()
    return H3_end + 1

def kink_anion(pose,abinfo):
    resi = abinfo.get_CDR_loop(h3).stop() - 1
    res = pose.residue(resi)
    atoms = []
    print "H3_N-1 (%i): %s" % (resi,res.name3())
    if res.name1() == "D":
        atoms.append(res.xyz("OD1"))
        atoms.append(res.xyz("OD2"))
    return atoms

def kink_cation(pose,abinfo):
    resi = abinfo.get_CDR_loop(h3).start() - 1
    res = pose.residue(resi)
    print "H3_0   (%i): %s" % (resi,res.name3())
    atoms = []
    if res.name1() == "R":
        atoms.append(res.xyz("NH1"))
        atoms.append(res.xyz("NH2"))
    if res.name1() == "K":
        atoms.append(res.xyz("NZ"))
    return atoms


def dihedral(p1,p2,p3,p4):
    a = p2-p1
    b = p3-p2
    c = p4-p3
    a.normalize()
    b.normalize()
    c.normalize()
    x = -a.dot(c) + a.dot(b) * b.dot(c)
    y = a.dot( b.cross(c) )
    angle = math.atan2(y,x)  # angle in [ -pi, pi ]  # degenerate if x or y == 0?
    angle = angle * 180 / math.pi
    return angle


def antibody_kink_Hbond(pose,abinfo):
    Aatoms = kink_anion(pose,abinfo)
    Catoms = kink_cation(pose,abinfo)

    HBdist = 100.0
    for Aa in Aatoms:
        for Ca in Catoms:
            HBdist = min( HBdist, (Aa - Ca).norm)
    if HBdist == 100.0: HBdist = 0

    return HBdist


def antibody_kink_bb_Hbond(pose,abinfo):

    Di = abinfo.get_CDR_loop(h3).stop() - 1
    D  = pose.residue(Di)
    DN = D.xyz("N")

    Ri = abinfo.get_CDR_loop(h3).start() - 1
    R  = pose.residue(Ri)
    RO = R.xyz("O")

    print "H3_DN   (%i): %s - %s" % (Di,D.name3(),DN)
    print "H3_RO   (%i): %s - %s" % (Ri,R.name3(),RO)

    bbHBdist = ( DN - RO ).norm

    return bbHBdist

def antibody_kink_trp_Hbond(pose,abinfo):

    Wi = abinfo.get_CDR_loop(h3).stop() + 1
    W  = pose.residue(Wi)
    if W.name3() != "TRP":
        return 0.0
    W_NE1 = W.xyz("NE1")

    kb1 = kink_begin(abinfo)
    kb  = pose.residue(kb1)
    kb_O = kb.xyz("O")

    WHBdist = ( W_NE1 - kb_O ).norm

    return WHBdist


def antibody_kink_geometry(pose, abinfo, debug=False):

    kb=kink_begin(abinfo)
    kr0=pose.residue(kb)
    kr1=pose.residue(kb+1)
    kr2=pose.residue(kb+2)
    kr3=pose.residue(kb+3)
    kseq = kr0.name1() + kr1.name1() + kr2.name1() + kr3.name1()

    print "Kink is defined from pose residues %i-%i: %s" % (kb,kb+3,kseq)

    if debug:
        pinfo = pose.pdb_info()
        print kr0, kr1, kr2, kr3
        print pinfo.number(kb),pinfo.number(kb+1),pinfo.number(kb+2),pinfo.number(kb+3)

    CA0=kr0.xyz("CA")
    CA1=kr1.xyz("CA")
    CA2=kr2.xyz("CA")
    CA3=kr3.xyz("CA")

    q = (CA0 - CA3).magnitude
    qbase = dihedral(CA0,CA1,CA2,CA3)

    return (q,qbase)


def main(args):
    '''Calculate kink geometry for a set of antibodies.  Usage: python kink_geom.py [pdb_files]
    '''
    if len(args) < 2:
        print
        print main.__doc__
        return

    outf = file('kink_geom.dat', 'w')
    outf.write("file       \t%10s\t%10s\t%10s\t%10s\t%10s\n" % ("q","qbase","HBdist","bbHBdist","W_HBdist") )
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

        q = antibody_kink_geometry(pose,abinfo)
        HBdist   = antibody_kink_Hbond(pose,abinfo)
        bbHBdist = antibody_kink_bb_Hbond(pose,abinfo)
        trpHBdist = antibody_kink_trp_Hbond(pose,abinfo)

        print "%s:\tq = %.3f, qbase = %.3f degrees, HBond_dist = %.3f Angstrom, bb_Hbond_dist = %.3f Angstrom, trp_HBond_dist = %.3f Angstrom" \
            % (filename,q[0],q[1],HBdist,bbHBdist,trpHBdist)
        outf.write( "%s\t%10.3f\t%10.3f\t%10.3f\t%10.3f\t%10.3f\n" % (filename,q[0],q[1],HBdist,bbHBdist,trpHBdist) )


if __name__ == "__main__": main(sys.argv)

