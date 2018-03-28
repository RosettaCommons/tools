from pdb_structure import *
from find_neighbors import *
from sequence_profile import *
import amino_acids as aas

def benchmark_dir( testset ) :
#    return "/home/andrew/andrew/rosetta/seqrec_bench/inputs/" + testset + "/"
    return "/Users/andrew/rosetta/seqrec_bench/" + testset + "/"

class SeqComp :
    def __init__( self ):
        self.rec_rate = 0.0
        self.npos = 0
        self.nrecd = 0
        self.test_pos_map = {}
        self.depth_classifier = DepthClassifier()
        self.aa_depth_recrates = []
        for d in range(self.depth_classifier.n_levels()) :
            self.aa_depth_recrates.append( {} )
            for aa in aas.amino_acids :
                self.aa_depth_recrates[ d ][ aa ] = [ 0, 0 ] # ordered pair: n recovered, n total
        self.total_aa_recrates = {}
        for aa in aas.amino_acids :
            self.total_aa_recrates[ aa ] = [ 0, 0 ] #ordered pair: n recovered, n total

    def add_position( self, ch, res, nneighbs, aa3, recovered ) :
        chres = ch + " " + res
        self.test_pos_map[ chres ] = ( nneighbs, recovered )
        if recovered :
            self.nrecd += 1.0
        self.npos += 1.0
        self.rec_rate = self.nrecd / self.npos
        if aa3 not in aas.longer_names : return
        aa1 = aas.longer_names[ aa3 ]
        d = self.depth_classifier.classify_depth_level( nneighbs )
        drecratelist = self.aa_depth_recrates[ d ][ aa1 ]
        recratelist = self.total_aa_recrates[ aa1 ]
        if recovered :
            drecratelist[0] += 1
            recratelist[0] += 1
        drecratelist[1] += 1
        recratelist[1] += 1

    def rec_rate_for_nneighb_range( self, nneighbs_low, nneighbs_high ) :
        count = 0
        ncorrect = 0
        for chres in list(self.test_pos_map.keys()):
            posdata = self.test_pos_map[ chres ]
            if posdata[0] <= nneighbs_high and posdata[0] >= nneighbs_low :
                count += 1
                if posdata[1]:
                    ncorrect += 1
        return ncorrect, count

    def compare_two_structures( self, ref_pdbstr, test_pdbstr, subset ) :
        # First determine the number of neighbors for the reference structure
        # which may take a little time.

        # reset
        self.rec_rate = 0
        self.npos = 0
        self.test_pos_map = {}

        ref_neighb_count = count_nneighbs_wi_cbeta_cutoff( ref_pdbstr, 10 );
        for ch,res in subset :
            self.add_position( ch, res, ref_neighb_count[ ch+" "+res ], ref_pdbstr.residue(ch,res).resname,
                               ref_pdbstr.residue( ch, res ).resname == test_pdbstr.residue( ch, res ).resname )

    # take a subset if you want, but if none is provided, then compare the entire structure
    def compare_two_pdbs( self, pdb1name, pdb2name, subset_file=None ):
        refpdb = PDBStructure()
        despdb = PDBStructure()
        refpdb.read_from_lines( open( pdb1name ).readlines() )
        despdb.read_from_lines( open( pdb2name ).readlines() )
        if subset_file :
            subset = load_subset( subset_file )
        else :
            subset = subset_from_pdb( refpdb )
        self.compare_two_structures( refpdb, despdb, subset )


def seqrec_for_subset( refpdb, despdb, subset ):
    count_total = 0
    same = 0
    for ch,res in subset :
        count_total += 1
        if refpdb.residue( ch, res ).resname == despdb.residue( ch, res ).resname :
            same += 1
    return same, count_total, float(same) / count_total

def compare_two_pdbs( pdb1name, pdb2name, subset_file=None ):
    refpdb = PDBStructure()
    despdb = PDBStructure()
    refpdb.read_from_lines( open( pdb1name ).readlines() )
    despdb.read_from_lines( open( pdb2name ).readlines() )
    if subset_file :
        subset = load_subset( subset_file )
    else :
        subset = subset_from_pdb( refpdb )
    return seqrec_for_subset( refpdb, despdb, subset )

def pairs_from_lines( lines ):
    pairs = []
    for line in lines :
        cols = line.split(" ")
        pairs.append( (cols[0], cols[1].strip() ))
    return pairs

def compare_all_pdbs( bench_set, verbose=True, swap_nats_and_designs=False ) :
    pdblist = [ x.strip() for x in open( benchmark_dir( bench_set ) + "pdb_names.list" ) ]
    recovery_list = []
    avg = 0
    count = 0

    for pdb in pdblist:
        nat_pdbname = benchmark_dir( bench_set ) + pdb + ".pdb"
        test_pdbname = pdb + ".pdb"

        if swap_nats_and_designs :
            # swap!
            temp = nat_pdbname
            nat_pdbname = test_pdbname
            test_pdbname = temp

        seq_comp = SeqComp()
        seq_comp.compare_two_pdbs( nat_pdbname, test_pdbname )
        recovery_list.append( (pdb, seq_comp) )

        if verbose:
            print(pdb, seq_comp.nrecd, seq_comp.npos, seq_comp.rec_rate)
            avg += seq_comp.rec_rate
            count += 1

    if verbose:
        print("average rec rate:",avg/count)
    return recovery_list
