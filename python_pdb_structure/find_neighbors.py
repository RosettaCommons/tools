from pdb_structure import *
#from cgkit.cgtypes import vec3
from vector3d import vector3d
from intdef_basics import *
import math

# Define an interface subset based on whether residues make inter-chain contacts
def find_interface_residues( pdb_structure, cutoff=10 ) :
    neighbs = find_neighbors_within_CB_cutoff( pdb_structure, cutoff )
    return find_interface_residues_from_neighbs( pdb_structure, neighbs )


def find_interface_residues_from_neighbs( pdb_structure, neighbs ) :
    intsets = {}
    for c in list(pdb_structure.chainmap.keys()):
        intsets[ c ] = set( [] )
    for neighb in neighbs:
        if neighb[0][0] != neighb[1][0] :
            intsets[ neighb[0][0] ].add( neighb[0][1] )
            intsets[ neighb[1][0] ].add( neighb[1][1] )
    return intsets

def cbeta_vector( pdb_structure, chain, resstring, n, ca, c, cb ):
    r = pdb_structure.residue( chain, resstring )
    if not r.has_atom( ca ) :
        #print "no ca:",ca, r.resstring, r.resname
        return None
    if r.has_atom( cb ) :
        cbv = (r.atom( cb ).xyz - r.atom( ca ).xyz)
        cbv.normalize()
        return cbv
    else :
        #print "no cb:",cb, r.resstring, r.resname
        if not r.has_atom( n ) or not r.has_atom( c )  :
            return None
        cac = (r.atom( ca ).xyz - r.atom( c ).xyz)
        can = (r.atom( ca ).xyz - r.atom( n ).xyz)
        cac.normalize()
        can.normalize()
        gly_cbv = cac + can;
        gly_cbv.normalize()
        return gly_cbv

def res1_pointed_at_res2( pdb_structure, c1, res1, c2, res2, n, ca, c, cb ):
    #print "res1 pointed at res2?",res1,res2
    cbv = cbeta_vector( pdb_structure, c1, res1, n, ca, c, cb )
    base = select_coord_for_residue( pdb_structure.residue( c1, res1 ), cb, ca)
    dest = select_coord_for_residue( pdb_structure.residue( c2, res2 ), cb, ca)
    #print cbv, base, dest
    if cbv and base and dest :
        dmb = (dest - base)
        dmb.normalize()
        costheta = cbv * dmb
        if costheta > math.cos( math.radians( 75 )) :
            return True
    return False

# return true if there are any sidechain atoms on at1 that are within
# a cutoff distance of any atoms on res2, or if res1 is gly and CA is
# within the cutoff distance of res2
def any_atoms_within_cutoff( res1, res2, cutoff ) :
    bbats = set( [ "CA", "N", "C", "O" ] )
    #for a1name in bbats:
    #    if res1.has_atom(a1name) :
    #        a1 = res1.atom(a1name)
    #        for a2name in bbats :
    #            if res2.has_atom(a2name) :
    #                a2 = res2.atom(a2name)
    #                print a1name,a2name,a1.xyz.x_,a1.xyz.y_,a1.xyz.z_,a2.xyz.x_,a2.xyz.y_,a2.xyz.z_
    #                if (a1.xyz - a2.xyz).length() < cutoff :
    #                    return True
    r1isgly = res1.resname == "GLY"
    for a1name in list(res1.stripped_atmap.keys()) :
        if not r1isgly and a1name in bbats: continue
        a1 = res1.stripped_atmap[a1name]
        for a2name in list(res2.stripped_atmap.keys()) :
            a2 = res2.stripped_atmap[a2name]
            if (a1.xyz - a2.xyz).length() < cutoff :
                return True
    return False

# compute the subset of residues that form disulfide bonds
def find_disulfided_residues_from_neighbs( pdb, neighbs ) :
    dslf_res = set( [] )
    for neighb in neighbs :
        c1,r1 = neighb[0]
        c2,r2 = neighb[1]
        res1, res2 = pdb.residue(c1,r1), pdb.residue(c2,r2)
        #print res1.resname, res2.resname
        if res1.resname == "CYS" and res2.resname == "CYS" :
            #print c1,r1,c2,r2, "are close"
            #n1,n2 = c1 + " " + r1, c2 + " " + r2
            if neighb[0] in dslf_res or neighb[1] in dslf_res :
                continue
            if not res1.has_atom( " SG " ) :
                continue
            if not res2.has_atom( " SG " ) :
                continue
            d = (res1.atom(" SG ").xyz - res2.atom( " SG " ).xyz).length()
            #print d
            if d < 2.4 :
                dslf_res.add( neighb[0] )
                dslf_res.add( neighb[1] )
    return dslf_res

def find_interface_pointing_residue( pdb ) :
    neighbs = find_neighbors_within_CB_cutoff( pdb, 12 )
    return find_interface_pointing_residues_from_neighbs( pdb, neighbs )

def find_interface_pointing_residues_from_neighbs( pdb, neighbs ):
    c  = " C  "
    ca = " CA "
    n  = " N  "
    cb = " CB "
    #cos70 = math.cos( math.radians( 70 ))
    intsets = {}
    for c in list(pdb.chainmap.keys()):
        intsets[ c ] = set( [] )
    for neighb in neighbs:
        if neighb[0][0] != neighb[1][0] :
            c1 = neighb[0][0]
            c2 = neighb[1][0]
            res1 = neighb[0][1]
            res2 = neighb[1][1]
            if res1 in intsets[ c1 ] and res2 in intsets[ c2 ] :
                continue
            #print "examining pair",c1,res1,c2,res2
            r1sc_near_r2 = any_atoms_within_cutoff( pdb.residue(c1,res1), pdb.residue(c2,res2), 5.5 )
            r2sc_near_r1 = any_atoms_within_cutoff( pdb.residue(c2,res2), pdb.residue(c1,res1), 5.5 )
            if r1sc_near_r2 :
                intsets[c1].add(res1)
                #print "adding res1"
            if r2sc_near_r1 :
                intsets[c2].add(res2)
                #print "adding res2"

            #print "cbeta neighbors?",res1,res2
            if res1 not in intsets[c1] :
                r1atr2 = res1_pointed_at_res2( pdb, c1, res1, c2, res2, n, ca, c, cb )
                if r1atr2 :
                    intsets[c1].add(res1)
                    #print "cbeta vector adding res1"
            if res1 not in intsets[c1] :
                r1atr2 = res1_pointed_at_res2( pdb, c2, res2, c1, res1, n, ca, c, cb )
                if r1atr2 :
                    intsets[c2].add(res2)
                    #print "cbeta vector adding res2"

            
    return intsets


def residue_center( r ) :
    center = vector3d()
    for at in r.atoms :
        center += at.xyz
    center.scale( 1.0 / len(r.atoms) )
    return center

def select_closest_atom( xyz1, residue ) :
    closest = None
    smallest_dist2 = -1
    for at in residue.atoms :
        d2 = at.xyz.distance_squared( xyz1 )
        if not closest or d2 < smallest_dist2 :
            closest = at
            smallest_dist2 = d2
    return at, smallest_dist2

def select_coord_for_residue( r, cb=" CB ", ca=" CA " ) :
    if r.has_atom( cb ) :
        coord = r.atom( cb ).xyz
    elif r.has_atom( ca ) :
        coord = r.atom( ca ).xyz
    else:
        # compute the mean coordinate of the residue
        coord = residue_center( r )
    return coord

# Return a list of residue pair structs (tuple of chainID + resstring)
# representing all residue pair cutoffs within a certain CBeta
# (CA for GLY) distance.  Naive O(N^2) implementation
def find_neighbors_within_CB_cutoff( pdb_struct, cutoff ) :
    cb = " CB "
    ca = " CA "
    neighbs = []
    for cind in range( len( pdb_struct.chains )):
        c = pdb_struct.chains[ cind ]
        for rind in range( len( c.residues )):
            r = c.residues[ rind ]
            coord = select_coord_for_residue( r, cb, ca )
            if not coord:
                continue
            for rind2 in range( rind+1, len(c.residues)):
                r2 = c.residues[ rind2 ]
                coord2 = select_coord_for_residue( r2, cb, ca )
                if not coord2:
                    continue
                #print cind, rind, cind, rind2, coord, coord2, (coord - coord2).length()
                if (coord - coord2).length() < cutoff :
                    neighbs.append( (( c.chain_name, r.resstring ), (c.chain_name, r2.resstring )) )
            for cind2 in range( cind + 1, len( pdb_struct.chains ) ):
                c2 = pdb_struct.chains[ cind2 ]
                for rind2 in range( len( c2.residues )):
                    r2 = c2.residues[ rind2 ]
                    coord2 = select_coord_for_residue( r2, cb, ca )
                    if not coord2:
                        continue
                    #print cind, rind, cind2, rind2, coord, coord2, (coord - coord2).length()
                    if (coord - coord2).length() < cutoff :
                        neighbs.append(( ( c.chain_name, r.resstring ), (c2.chain_name, r2.resstring ) ))
    return neighbs

def count_nneighbs_wi_cbeta_cutoff( pdb_struct, cutoff, neighbs=None ):
    if neighbs==None:
        neighbs = find_neighbors_within_CB_cutoff( pdb_struct, cutoff )
    neighb_count = {}
    for neighb in neighbs :
        ch1res1 = neighb[0][0]+" "+neighb[0][1]
        ch2res2 = neighb[1][0]+" "+neighb[1][1]
        if ch1res1 not in neighb_count :
            neighb_count[ch1res1] = 0
        if ch2res2 not in neighb_count :
            neighb_count[ch2res2] = 0
        neighb_count[ch1res1] += 1
        neighb_count[ch2res2] += 1
    return neighb_count

def neighbor_graph( pdb_struct, cutoff=10, neighbs=None ) :
    if neighbs==None:
        neighbs = find_neighbors_within_CB_cutoff( pdb_struct, 10 )
    neighb_graph = {}
    for neighb in neighbs :
        c1,r1=neighb[0]
        c2,r2=neighb[1]
        if c1 not in neighb_graph : neighb_graph[ c1 ] = {}
        if c2 not in neighb_graph : neighb_graph[ c2 ] = {}
        if r1 not in neighb_graph[ c1 ] : neighb_graph[c1][r1] = []
        if r2 not in neighb_graph[ c2 ] : neighb_graph[c2][r2] = []
        neighb_graph[c1][r1].append( (c2,r2) )
        neighb_graph[c2][r2].append( (c1,r1) )
    return neighb_graph

def get_average_neighb_vector( pdb_struct, target_res, neighbor_cutoff=10, close_contact_max_dist=10, target_chain="all", neighb_graph=None ):
    avg_vect = vector3d()
    if neighb_graph==None :
        neighbs = find_neighbors_within_CB_cutoff( pdb_struct, 10 )
        neighb_graph = neighbor_graph( pdb_struct, 10, neighbs )
    if target_chain == "all" :
        chains = list(pdb_struct.chainmap.keys())
    else :
        chains = [ target_chain ]
    target_xyz = select_coord_for_residue( target_res )
    ccd2 = close_contact_max_dist * close_contact_max_dist;
    for ch in chains :
        for neighb in neighb_graph[ target_res.chain.chain_name ][ target_res.resstring ] :
            if neighb[0] == ch :
                closest_neighb, d2 = select_closest_atom( target_xyz, pdb_struct.residue( neighb[0], neighb[1] ) )
                if d2 > ccd2 : continue #skip residues that are too far away
                vect = closest_neighb.xyz - target_xyz
                vect.normalize()
                avg_vect += vect
    avg_vect.normalize()
    return avg_vect
