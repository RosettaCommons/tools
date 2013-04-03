#!/usr/bin/python

from sys import stdout, stderr
import string
from math import sqrt
from parse_options import parse_options

def get_positions( lines, atom_names_in ):

    count = 0
    positions = []
    which_res = []
    atom_names = []
    oldresnum = ''
    for line in lines:

        if (len(line)>54 and  line[0:4] == 'ATOM' ):

            resnum = line[23:26]
            if not resnum == oldresnum:
                count = count + 1

            atom_name = line[12:16]
            position = [float(line[30:38]),float(line[38:46]),float(line[46:54])]

            if atom_name in atom_names_in:
                positions.append( position )
                which_res.append( count )
                atom_names.append( atom_name )

            oldresnum = resnum

    return ( positions, which_res, atom_names, count )

def get_dist( pos1, pos2 ):
    dist2 = 0.0
    for k in range(3):
        dsub = ( pos1[k] - pos2[k] )
        dist2 += dsub*dsub

    return sqrt( dist2 )

def write_function( fid, dist, fade, stdev_to_use, cst_depth, tight_fade ):
    if fade:
        stdout.write( " FADE" )
        if tight_fade:
            stdout.write( " %8.3f %8.3f %8.3f" %
                      ( dist - 2 * stdev_to_use, dist + 2 * stdev_to_use, 2*stdev_to_use ) )
        else:
            stdout.write( " %8.3f %8.3f %8.3f" %
                           ( dist - 2*stdev_to_use, dist + 2*stdev_to_use, stdev_to_use ) )

        if not ( cst_depth == 0.0 ):
            stdout.write( " %8.2f 0.0 \n" % cst_depth )
        else:
            stdout.write( " -10.0 10.0 \n" )

    else:
        assert( not tight_fade )
        if( cst_depth != 0.0 ):
            print 'should not supply cst_depth with harmonic constraints!'
            exit( 0 )
        stdout.write( " HARMONIC %8.3f %8.3f \n" % \
                          ( dist,stdev_to_use) )


def generate_constraints( argv, atom_names_in1, atom_names_in2, dist_cut_default ):

    pdbfile = argv[1]
    fade = parse_options( argv, "fade", 0 )
    DIST_CUT = parse_options( argv, "dist_cut", dist_cut_default )
    STDEV = parse_options( argv, "stdev", 2.0 )
    STDEV_LOOSE = parse_options( argv, "stdev_loose", 4.0 )
    SEQ_SEP_CUTOFF = parse_options( argv, "seq_sep_cutoff", 3 )
    loose_res = parse_options( argv, "loose_res", [-1] )
    cst_res = parse_options( argv, "cst_res", [-1] )
    no_cst_res = parse_options( argv, "no_cst_res", [-1] )
    COORD_CST = parse_options( argv, "coord_cst", 0 )
    anchor_atom_name = parse_options( argv, "anchor_atom", " CA " )
    anchor_resnum = parse_options( argv, 'anchor_res', 0 )
    cst_depth = parse_options( argv, 'cst_depth', 999.9 )
    tight_fade = parse_options( argv, 'tight_fade', 0 )
    if (cst_depth == 999.9) : cst_depth = 0.0

    if ( DIST_CUT == 0.0 and len( atom_names_in2) == 0 ): COORD_CST = 1 #signal to use coordinate constraints

    if ( COORD_CST and anchor_resnum == 0 ):
        stderr.write( 'If you want coord_cst, need to also specify -anchor_res\n' )
        exit( 0 )

    if ( not isinstance( atom_names_in1, list ) ): atom_names_in1 = [ atom_names_in1 ]
    if ( not isinstance( atom_names_in2, list ) ): atom_names_in2 = [ atom_names_in2 ]

    lines = open( pdbfile ).readlines()
    ( positions1, which_res1, atom_names1, totres) = get_positions( lines, atom_names_in1 )
    ( positions2, which_res2, atom_names2, totres) = get_positions( lines, atom_names_in2 )

    cst_file = pdbfile+'.cst'

    fid = open( 'TEST.pml', 'w' )
    fid.write( 'reinitialize\n' )
    fid.write( 'load %s, model_pdb\n' % pdbfile )
    fid.write('hide everything, model_pdb\n' )
    fid.write('show cartoon, model_pdb\n' )
    fid.write('cmd.spectrum(selection = "model_pdb")\n' )
    fid.write('bg_color white\n' )

    if len( cst_res ) == 0:
        for m in range( 1, totres+1):
            if m not in no_cst_res:
                cst_res.append( m )

    if len( loose_res ) > 0:
        command = 'color gray, resi %d' % loose_res[0]
        for m in loose_res[1:]: command += '+%d' % m
        fid.write( command + '\n' )

    free_res = []
    for i in range( 1, totres+1 ):
        if ( i not in cst_res ) and ( i not in loose_res ):
            free_res.append( i )

    if len( free_res ) > 0:
        command = 'color white, resi %d' % free_res[0]
        for m in free_res[1:]: command += '+%d' % m
        fid.write( command + '\n' )

    if COORD_CST:
        stdout.write( "[ coordinates ]\n" )

        for i in range( len( which_res1) ):

            res1 = which_res1[ i ]
            if ( res1 in free_res ): continue
            atom1 = atom_names1[ i ]

            stdev_to_use = STDEV
            if res1 in loose_res:  stdev_to_use = STDEV_LOOSE

            stdout.write( " %s %d   %s %d    %8.3f %8.3f %8.3f" % \
                           ( atom1, res1,
                             anchor_atom_name,
                             anchor_resnum,
                             positions1[ i ][0],
                             positions1[ i ][1],
                             positions1[ i ][2] ) )

            write_function( stdout, 0.0, fade, stdev_to_use, cst_depth, tight_fade )

    else: # AtomPairs.
        stdout.write( "[ atompairs ]\n" )
        connect_pairs = []
        for i in range( len( which_res1) ):

            res1 = which_res1[ i ]
            if ( res1 in free_res ): continue
            atom1 = atom_names1[ i ]

            for j in range( len( which_res2) ):

                res2 = which_res2[ j ]
                if ( res2 in free_res ): continue
                atom2 = atom_names2[ j ]

                if ( abs( res1 - res2 ) <= SEQ_SEP_CUTOFF ): continue

                dist = get_dist( positions1[ i ], positions2[ j ] )

                if ( dist < DIST_CUT ):

                    if ( [ [res2, atom2], [res1, atom1] ] in connect_pairs ): continue #no redundancy

                    connect_pairs.append( [ [ res1, atom1 ], [res2, atom2 ] ] )

                    if ( ( res1 in loose_res ) or ( res2 in loose_res ) ):
                        stdev_to_use = STDEV_LOOSE
                        fid.write('dist LOOSE_CST, resi %d and name %s, resi %d and name %s\n' % (res1, atom1, res2, atom2 ))
                    else:
                        stdev_to_use = STDEV;
                        fid.write('dist CST, resi %d and name %s, resi %d and name %s\n' % (res1, atom1, res2, atom2 ))

                    stdout.write( " %s %d   %s %d  " %   ( atom1, res1, atom2, res2 ) )
                    write_function( stdout, dist, fade, stdev_to_use, cst_depth, tight_fade )

        fid.write( 'hide labels, CST\n')
        fid.write( 'hide labels, LOOSE_CST\n')
        fid.write( 'color red, CST\n')
        fid.write( 'color gray, LOOSE_CST\n')

    if len( loose_res ) > 0 :
        fid.write( 'color gray, resi %d' % loose_res[0] )
        for m in loose_res[1:]: fid.write('+%d'%m)
        fid.write('\n')

    fid.close()
    stderr.write( 'Made a pymol script in TEST.pml\n' )
