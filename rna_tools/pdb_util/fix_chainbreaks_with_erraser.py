#!/usr/bin/python

from sys import argv
from os import popen,system
from parse_options import get_ints
import string

pdbs = argv[1:]

for pdb in pdbs:
    resnum_line = popen( 'get_res_num.py  '+pdb ).readlines()[0][:-1]
    cols = resnum_line.split()
    resnum = []
    print resnum_line
    for col in cols:   get_ints( col, resnum )
    print 'Residue numbers: ', resnum

    # find existing cutpoint
    command = 'rna_rosetta_ready_set.py -pdb %s' % pdb
    lines = popen( command ).readlines()
    cutpoint_num = []
    for line in lines:
        if ( line.find( 'cutpoint in Rosetta pdb file' ) > -1 ):
            cols = line[:-1].split( ':' )[-1].split()
            for col in cols: cutpoint_num.append( int( col) )

    print 'Cutpoint number in renumbered PDB according to rna_rosetta_ready_set:', cutpoint_num

    # actually use my new script 'check_cutpoint.py' to figure out stuff.
    command = 'check_cutpoints.py ' + pdb
    lines = popen( command ).readlines()
    cutpoints_num        = map( lambda x:int(x), string.split( string.split( lines[1][:-1], ':'  )[1]) )
    distort_geometry_num = map( lambda x:int(x), string.split( string.split( lines[2][:-1], ':'  )[1]) )

    print 'Cutpoint number in renumbered PDB according to check_cutpoint.py:', cutpoint_num
    print 'Distorted C5*-O5* geometry in renumbered PDB according to check_cutpoint.py:', distort_geometry_num

    rebuild_num  = distort_geometry_num
    for n in cutpoint_num:
        if ( resnum[n] > resnum[n-1]+1 ) :
            print n, ' is meant to be a cutpoint, so no rebuild there.'
            continue # its meant to be a cutpoint, actually
        if ( n   not in rebuild_num ): rebuild_num.append( n )
        if ( n+1 not in rebuild_num ): rebuild_num.append( n+1 )

    print rebuild_num

    tag = pdb.replace( '.pdb', '')

    if len( rebuild_num ) > 0:
        print rebuild_num
        rebuild_num_list = ''
        for m in rebuild_num: rebuild_num_list += ' %d' % m

        readysetpdb = tag+'_ready_set.pdb'

        command = 'seq_rebuild.py -pdb  %s -scoring_file rna/rna_hires_07232011_with_intra_base_phosphate  -rebuild_res_list  %s  -native_screen_RMSD 2.0' % (readysetpdb, rebuild_num_list)
        print command
        system( command )

        command = 'cp  %s_ready_set_erraser.pdb   %s_fixcutpoint.pdb' % (tag,tag)
        print command
        system( command )

        command = 'renumber_pdb_in_place.py %s_fixcutpoint.pdb   %s' (tag,resnum_line)
        print command
        system( command )
    else:
        print 'No cutpoints identified!!!!'
        command = 'cp %s.pdb  %s_fixcutpoint.pdb' % (tag,tag)
        print command
        system( command )

