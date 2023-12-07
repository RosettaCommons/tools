import pdb_structure
import amino_acids

def print_alignments( a1, a2, name1, name2 ) :
    # wrap at 60 chars
    assert( len(a1) == len(a2) )
    cnt1, cnt2 = 0, 0
    seqstr = [None]*60;
    while cnt1 < len(a1) :
        for i in range( 60 ) :
            seqstr[i] = a1[cnt1]
            cnt1 += 1
            if cnt1 == len(a1) : 
                seqstr = seqstr[:i+1]
                break
        print(name1, "".join(seqstr))
        for i in range( 60 ) :
            seqstr[i] = a2[cnt2]
            cnt2 += 1
            if cnt2 == len(a2) : 
                break
        print(name2, "".join(seqstr))
        print()


def print_mat( mat ) :
    for i in range(len(mat)) :
        print("[", end=' ')
        for j in range(len(mat[i])) :
            print("%3d" % mat[i][j], end=' ')
        print("]")

def needleman_wunsch( seq1, seq2 ) :
    gap_pen = -1

    score = [ None ] * ( len(seq1)+1 )
    trace = [ None ] * ( len(seq1)+1 )

    # three moves:
    # 1: increase i, keep j the same  = insert a gap into seq2
    # 2: increase j, keep i the same  = insert a gap into seq1
    # 3: increase i and j             = align i and j to each other

    for i in range(len(seq1)+1) :
        score[i] = [ None ] * ( len(seq2) +1 )
        trace[i] = [ None ] * ( len(seq2) +1 )
        score[i][0] = gap_pen * i
        trace[i][0] = 1
    for j in range(len(seq2)+1) :
        score[0][j] = gap_pen * j
        trace[0][j] = 2

    for i in range(1,len(seq1)+1) :
        for j in range(1,len(seq2)+1) :
            sc1 = gap_pen + score[i-1][j]
            sc2 = gap_pen + score[i][j-1]
            sc3 = ( 1 if seq1[i-1] == seq2[j-1] else 0 ) + score[i-1][j-1]
            which = 1 if ( sc1 >= sc2 and sc1 >= sc3 ) else ( 2 if ( sc2 >= sc1 and sc2 >= sc3 ) else 3 )
            # print i, j, seq1[i-1], seq2[j-1], sc1, sc2, sc3, which
            trace[i][j] = which
            score[i][j] = sc1 if which == 1 else ( sc2 if which == 2 else sc3 )
            #print "i",i,"j", j,"si", seq1[i-1], "sj",seq2[j-1], "sc1", sc1, "sc2",sc2,"sc3", sc3, "trij", trace[i][j], "sci-1j", score[i-1][j], "scij-1", score[i][j-1], "sci-1j-1", score[i-1][j-1],"scij", score[i][j]

    #print "score"
    #print_mat( score )
    #print "trace"
    #print_mat( trace )

    # backtrace
    best_path = []
    i = len(seq1)
    j = len(seq2)
    while i != 0 and j != 0 :
        #print "best path, appending ", i, j, "with", trace[i][j]
        best_path.append( trace[i][j] )
        if trace[i][j] == 3 :
            i, j = i-1, j-1
        elif trace[i][j] == 2 :
            i, j = i, j-1
        elif trace[i][j] == 1 :
            i, j = i-1, j

    #print "at conclusion of backtrace", i, j

    while j != 0 :
        best_path.append( 2 )
        j = j-1
    while i != 0 :
        best_path.append( 1 )
        i = i-1


    best_path = best_path[::-1]

    seq1_aln, seq2_aln = [], []
    cnt1, cnt2 = 0, 0
    for i in range(len(best_path)) :
        if best_path[i] == 2 :
            seq1_aln.append( "-" )
            seq2_aln.append( seq2[cnt2] )
            cnt2 += 1
        elif best_path[i] == 1 :
            seq1_aln.append( seq1[cnt1] )
            seq2_aln.append( "-" )
            cnt1 += 1
        elif best_path[i] == 3 :
            seq1_aln.append( seq1[cnt1] )
            seq2_aln.append( seq2[cnt2] )
            cnt1 += 1
            cnt2 += 1
    return seq1_aln, seq2_aln

def oneletter_seq_for_chain( chain ) :
    seq = []
    for res in chain.residues : 
        if res.resname in amino_acids.longer_names :
            seq.append( amino_acids.longer_names[ res.resname ])
    return seq
