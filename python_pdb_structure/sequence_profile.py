import amino_acids as aas
import math

class SequenceProfile :
    def __init__( self ) :
        self.aa_counts = {}
        for aa in aas.amino_acids :
            self.aa_counts[ aa ] = 0
        self.tot_aa = 0
    def add_obs( self, aa ) :
        assert( aa in self.aa_counts )
        self.aa_counts[ aa ] += 1
        self.tot_aa += 1
    def aa_freq( self, aa ) :
        return float( self.aa_counts[ aa ] ) / float( self.tot_aa )
    def increment( self, sqprf2 ) :
        for aa in aas.amino_acids :
            self.aa_counts[ aa ] += sqprf2.aa_counts[ aa ]
        self.tot_aa += sqprf2.tot_aa

# default depth classification [0-16] neighbors exposed, [17-23] intermediate, [24+] buried
class DepthClassifier :
    def __init__( self, depth_names=None, upper_boundaries=None ) :
        if depth_names and upper_boundaries :
            assert( len(depth_names) == len(upper_boundaries)+1 )
            self.depth_names = depth_names
            self.upper_boundaries = upper_boundaries
        else :
            self.depth_names = [ "EXP", "INT", "BUR" ]
            self.upper_boundaries = [ 16, 23 ]
        self.max_depth_name_length = 3
        for name in self.depth_names :
            if len(name) > self.max_depth_name_length : self.max_depth_name_length = len(name)
    def classify_depth_level( self, n_neighbors ) :
        for i in range(len(self.upper_boundaries)):
            if n_neighbors <= self.upper_boundaries[i] :
                return i #found it
        return len(self.upper_boundaries) # return the index of the most buried depth
    def n_levels( self ) :
        return len(self.depth_names)

class DepthSequenceProfile :
    def __init__( self ) : # default ["EXP", "INT", "BUR" ]
        self.depth_profiles = []
        self.depth_classifier = DepthClassifier()
        for i in range( self.depth_classifier.n_levels ) :
            self.depth_profiles.append( SequenceProfile() )
        self.whole_profile = SequenceProfile()
    def classify_depth_level( self, n_neighbors ) :
        return self.depth_classifier( n_neighbors )
    def add_obs( self, aa, n_neighbors ) :
        self.whole_profile.add_obs( aa )
        depth_level = self.classify_depth_level( n_neighbors )
        #print "dsp: add_obs:", aa, n_neighbors, depth_level
        self.depth_profiles[ depth_level ].add_obs( aa )
        self.whole_profile.add_obs(aa)
    def write_output( self ) :
        strlines = []
        currline = " " * ( self.depth_classifier.max_depth_name_length + 1 )
        for aa in aas.amino_acids :
            currline += "     " + aa
        strlines.append( currline + "\n" )
        currline = "ALL" + " " * ( self.depth_classifier.max_depth_name_length + 1 - 3 )
        strlines.append( currline + self.strline_for_seqprof( self.whole_profile ) + "\n" )
        for i in range( len( self.depth_classifier.depth_names )) :
            strlines.append( self.depth_classifier.depth_names[ i ]
                             + " " * ( self.depth_classifier.max_depth_name_length - len( self.depth_classifier.depth_names[ i ] ) + 1 )
                             + self.strline_for_seqprof( self.depth_profiles[ i ] ) + "\n" )
        return strlines
    def strline_for_seqprof( self, seqprof ) :
        line = ""
        for aa in aas.amino_acids :
            line += " %5.1f" % ( 100 * seqprof.aa_freq( aa ) )
        return line
    def increment( self, sqprf2 ) :
        assert( len(self.depth_profiles) == len( sqprf2.depth_profiles))
        for i in range(len(self.depth_profiles)):
            self.depth_profiles[i].increment( sqprf2.depth_profiles[i] )
        self.whole_profile.increment( sqprf2.whole_profile )

# prof1 is considered the 'true' distribution, prof2 is the 'approximation'
def cross_entropy_for_twoseqprofs( prof1, prof2 ) :
    entropy = 0.0
    cross_entropy = 0.0
    for aa in aas.amino_acids :
        entropy -= prof1.aa_freq(aa) * math.log( prof1.aa_freq(aa) )
        log_prof2 = math.log( 0.00001 )
        if prof2.aa_freq(aa) != 0.0 :
            log_prof2 = math.log( prof2.aa_freq(aa) )
        cross_entropy -= prof1.aa_freq( aa ) * log_prof2
    return cross_entropy, entropy

def print_depth_cross_entropies( dprof1, dprof2 ) :
    xent, ent = cross_entropy_for_twoseqprofs( dprof1.whole_profile, dprof2.whole_profile )
    print("ALL :",( "%6.3f %6.3f %6.3f" % (ent, xent, (xent - ent) )))
    for i in range(len(dprof1.depth_names)) :
        xent, ent = cross_entropy_for_twoseqprofs( dprof1.depth_profiles[i], dprof2.depth_profiles[i] )
        print(dprof1.depth_names[i], ":", ( "%6.3f %6.3f %6.3f" % (ent, xent, (xent - ent) )))

def freqs_for_aa_group( seqprof, aaset ) :
    sum = 0
    for aa in aas.amino_acids :
        if aa in aaset :
            sum += seqprof.aa_counts[ aa ]
    return float( sum ) / seqprof.tot_aa
