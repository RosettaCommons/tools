from sequence_profile import *
from pdb_structure import *
from find_neighbors import *
import sys
import amino_acids as aa

#def depth_given_ncontacts( depths, ncontacts ) :
#    for i in xrange( len( depths ) ) :
#        if ncontacts >= depths[i][0] and ncontacts <= depths[i][1 ] :
#            return i
#    return len(depths)-1

def compute_depth_seqprof_for_structure( seqprofs, pdb, ncontacts ) :
    for ch in pdb.chains :
        for res in ch.residues :
            if res.resname in aa.longer_names:
                seqprofs.add_obs( aa.longer_names[ res.resname ], ncontacts[ res.resid() ] )


def compute_seqprof_for_list( listfile ) :
    seqprofs = DepthSequenceProfile()
    for pdbfile in listfile :
        if pdbfile[-3:] != ".pdb" :
            pdbfile += ".pdb"
            pdb = pdbstructure_from_file( pdbfile )
        ncontacts = count_nneighbs_wi_cbeta_cutoff( pdb, 10 )
        compute_depth_seqprof_for_structure( seqprofs, pdb, ncontacts )
    return seqprofs

def add_pdb_ext( pdbfile ) :
    if pdbfile[-3:] != ".pdb" :
        pdbfile += ".pdb"
    return pdbfile


def compute_seqprof_for_two_lists( listfile1, listfile2, depth_level_upper_limits=None ) :
    assert( len( listfile1 ) == len( listfile2 ) )
    seqprofs1 = DepthSequenceProfile()
    seqprofs2 = DepthSequenceProfile()
    if depth_level_upper_limits : # if not specified, use the defaults
        seqprofs1.upper_boundaries = depth_level_upper_limits
        seqprofs2.upper_boundaries = depth_level_upper_limits
    for i in range( len( listfile1 ) ):
        pdbfile1 = add_pdb_ext( listfile1[ i ] )
        pdbfile2 = add_pdb_ext( listfile2[ i ] )
        pdb1 = pdbstructure_from_file( pdbfile1 )
        pdb2 = pdbstructure_from_file( pdbfile2 )
        ncontacts = count_nneighbs_wi_cbeta_cutoff( pdb1, 10 )
        compute_depth_seqprof_for_structure( seqprofs1, pdb1, ncontacts )
        compute_depth_seqprof_for_structure( seqprofs2, pdb2, ncontacts )
    return seqprofs1, seqprofs2

if __name__ == "__main__" :
    #ranges = [ (0,16), (17,23), (24,100000) ]
    depth_level_upper_limits = [ 16, 23 ]
    if len( sys.argv ) == 3 :
        listfile1 = [ x.strip() for x in open( sys.argv[1] ).readlines() ]
        listfile2 = [ x.strip() for x in open( sys.argv[2] ).readlines() ]
        seqprofs1, seqprofs2 = compute_seqprof_for_two_lists( listfile1, listfile2, depth_level_upper_limits )
        seqproflines1 = seqprofs1.write_output()
        seqproflines2 = seqprofs2.write_output()
        print("Profile 1:")
        for line in seqproflines1 :
            print(line, end=' ')
        print("Profile 2:")
        for line in seqproflines2 :
            print(line, end=' ')
        print_depth_cross_entropies( seqprofs1, seqprofs2 )
        
