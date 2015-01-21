#!/usr/bin/python

from rosetta import *
init()

from rosetta.protocols.antibody import *
from numpy import *
import fnmatch
import os

def compute_pairwise_flex(allPoses,allAbInfos):
    npdbs = len(allPoses)
    flex = zeros(shape=(npdbs,npdbs))
    for p1,ai1,i in zip(allPoses,allAbInfos,range(npdbs)):
        print i,p1.pdb_info().name()
        for p2,ai2,j in zip(allPoses,allAbInfos,range(npdbs)):
            if j<i:
                flex[i,j] = 0
            else:
                flex[i,j] = global_loop_rmsd(p1,p2,ai2.get_CDR_in_loopsop(h3))

    print flex
    iu = triu_indices(npdbs, 1)
    upper_diagonal_mean = mean( flex[iu] )
    return upper_diagonal_mean


def load_pdbs(pdblist):
    poses = []
    abinfos = []
    for filename in pdblist:
        print filename

        try:
            pose = pose_from_pdb(filename)
            abinfo = AntibodyInfo(pose)
        except:
            print 'Problems reading %s, skipping.' % filename
            continue

        poses.append(pose)
        abinfos.append(abinfo)

    if len(poses) < 2:
        print "Loaded fewer than two poses...no flex to compute.  Exiting."
        return

    print poses[0]
    print abinfos[0]

    return poses, abinfos


def align_pdb_H_chains(poses,abinfos):
    # align pdbs to pose 0
    for pose, abinfo in zip(poses[1:],abinfos[1:]):
        align_to_native(pose,poses[0],abinfo,abinfos[0],"H") # pose is an OP so the align sticks


def main_pdblist(args):
    '''Usage: python H3flex.py [pdb_files]
       Calculate rmsd (flex) between a set of antibody H3 loops.
       Flex is defined as the average rmsd between all pairs of H3 loops in the set
       First pdb file will be used for alignment
    '''
    if len(args) < 2:
        print
        print main_pdblist.__doc__
        return

    pdblist = args[1:]
    poses,abinfos = load_pdbs(pdblist)
    align_pdb_H_chains(poses,abinfos)
    flex = compute_pairwise_flex(poses,abinfos)
    print 'Mean flex: ', flex
    return


def main_dirlist(args):
    '''Usage: python H3flex.py [directories]
       Calculate flex between a set of antibody H3 loops for a set of directories
       Flex is defined as the average rmsd between all pairs of H3 loops in the set
       Assumes pdb files are dir/top10/*pdb
       First pdb file will be used for alignment
    '''
    if len(args) < 1:
        print
        print main_dirlist.__doc__
        return

    dirlist = args[1:]
    startdir = os.getcwd()

    outfname = os.path.basename(startdir) + ".H3flex"
    outf = file(outfname, 'w')
    outf.write( "%s\t %10s %10s %10s\n" % ("Ab","H3len","Npdbs","H3flex") )
    for dir in dirlist:
        try:
            os.chdir(dir+"/top10")
            pdblist = fnmatch.filter( os.listdir("./"),"*pdb*" )
            print dir,":", pdblist

            poses,abinfos = load_pdbs(pdblist)
            align_pdb_H_chains(poses,abinfos)
            flex = compute_pairwise_flex(poses,abinfos)
            os.chdir(startdir)

            print dir,": ",flex
            outf.write( "%s\t %10i %10i %10.3f\n" % (dir,abinfos[0].get_CDR_length(h3),len(poses),flex) )
            outf.flush()
        except:
            print "Error processing directory ",dir
            os.chdir(startdir)
            continue
    return


if __name__ == "__main__": main_dirlist(sys.argv)

