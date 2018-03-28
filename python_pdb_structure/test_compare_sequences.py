from compare_sequences import *
from pdb_structure import *
import sys
import amino_acids as aas

# $1 == bench set

if len(sys.argv) < 2 :
    print("ERROR: benchmark set directory name must be specified.")
    sys.exit(1)

bench_set = sys.argv[1]
seq_comps = compare_all_pdbs( bench_set, True )
rev_comps = compare_all_pdbs( bench_set, True, True )

tot_rec, tot_n = 0.0, 0.0
tot_rec_exp, tot_n_exp, tot_rec_int, tot_n_int, tot_rec_bur, tot_n_bur = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
for pdb,seqrec in seq_comps :
    id_tot_rec, id_n_tot = 0,0
    id_rec_exp, id_n_exp, id_rec_int, id_n_int, id_rec_bur, id_n_bur = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    #for seqrec in ss[ intdef ]:

    id_tot_rec = seqrec.nrecd
    id_n_tot   = seqrec.npos
    rec_exp, n_exp = seqrec.rec_rate_for_nneighb_range(  0, 16 )
    rec_int, n_int = seqrec.rec_rate_for_nneighb_range( 17, 23 )
    rec_bur, n_bur = seqrec.rec_rate_for_nneighb_range( 24, 100 )
    id_rec_exp += rec_exp
    id_n_exp += n_exp
    id_rec_int += rec_int
    id_n_int += n_int
    id_rec_bur += rec_bur
    id_n_bur += n_bur

    n_ss_intdef = 1 #len( ss[ intdef ] )
    id_tot_rec /= n_ss_intdef
    id_n_tot /= n_ss_intdef
    id_rec_exp /= n_ss_intdef
    id_rec_int /= n_ss_intdef
    id_rec_bur /= n_ss_intdef
    id_n_exp /= n_ss_intdef
    id_n_int /= n_ss_intdef
    id_n_bur /= n_ss_intdef

    tot_rec += id_tot_rec
    tot_n += id_n_tot
    tot_rec_exp += id_rec_exp
    tot_rec_int += id_rec_int
    tot_rec_bur += id_rec_bur
    tot_n_exp += id_n_exp
    tot_n_int += id_n_int
    tot_n_bur += id_n_bur

print("Total recovered", tot_rec/tot_n, "(", tot_rec, "/", tot_n,  ")")
print("Exposed recovery", tot_rec_exp/tot_n_exp, "(", tot_rec_exp , "/", tot_n_exp ,  ")")
print("Intermediate recovered", tot_rec_int/tot_n_int, "(", tot_rec_int , "/", tot_n_int ,  ")")
print("Buried recovered", tot_rec_bur/tot_n_bur, "(", tot_rec_bur , "/", tot_n_bur ,  ")")

first_pdb,first_seqcomp = seq_comps[0]
rev_first_pdb,rev_first_seqcomp = rev_comps[0]
print("AA ", end=' ')
for name in first_seqcomp.depth_classifier.depth_names :
    print(name + " " * (23 -len(name)), end=' ')
print("TOT")
for aa in aas.amino_acids :
    print(aas.one_letter_names[aa], end=' ')
    for d in range(first_seqcomp.depth_classifier.n_levels() ) :
        count_rec = 0
        count_tot = 0
        for pdb,seqrec in seq_comps :
            count_rec += seqrec.aa_depth_recrates[d][aa][0]
            count_tot += seqrec.aa_depth_recrates[d][aa][1]
        print("( %5.3f, %5d, %5d )" % ( float(count_rec)/float(count_tot) , count_rec, count_tot ), end=' ')
    count_rec = 0
    count_tot = 0
    for pdb,seqrec in seq_comps :
        count_rec += seqrec.total_aa_recrates[aa][0]
        count_tot += seqrec.total_aa_recrates[aa][1]
    print("( %5.3f, %5d, %5d )" % ( float(count_rec)/float(count_tot) , count_rec, count_tot ))

print()
print("Design composition")
for aa in aas.amino_acids :
    print(aas.one_letter_names[ aa ], end=' ')
    for d in range( rev_first_seqcomp.depth_classifier.n_levels() ) :
        count_rec = 0
        count_tot = 0
        for pdb,seqrec in rev_comps :
            count_rec += seqrec.aa_depth_recrates[d][aa][0]
            count_tot += seqrec.aa_depth_recrates[d][aa][1]
        print("( %5.3f, %5d, %5d )" % ( float(count_rec)/float(count_tot), count_rec, count_tot ), end=' ')
    count_rec = 0
    count_tot = 0
    for pdb,seqrec in rev_comps :
        count_rec += seqrec.total_aa_recrates[aa][0]
        count_tot += seqrec.total_aa_recrates[aa][1]
    print("( %5.3f, %5d, %5d )" % ( float(count_rec)/float(count_tot) , count_rec, count_tot ))
