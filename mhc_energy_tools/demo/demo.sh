# ======================================== 
# Scoring

#1 List Propred epitope hits (5%) for SakStar, using the default Southwood set of 8 alleles

../score.py --fa in/2sak_A.fasta --report hits 

#2 Using native numbering, for ref against Choi Fig. 6

../score.py --fa in/2sak_A_native.fasta --report hits

#3 While the printout is useful from the command-line, the csv file is better for posterity and potentially for further analysis

../score.py --fa in/2sak_A.fasta --report hits --report_file out/2sak_A.pp5.csv

#4 Compare to NetMHCII (10% -- "weak binder" as of NetMHCII 2.3 web server), for the Greenbaum set of alleles

../score.py --fa in/2sak_A.fasta --report hits --netmhcii --epi_thresh 10 --allele_set greenbaum11 --report_file out/2sak_A.nm10g.csv

#5 And the Paul set

../score.py --fa in/2sak_A.fasta --report hits --netmhcii --epi_thresh 10 --allele_set paul15 --report_file out/2sak_A.nm10p.csv

#6 Theese things can be easier to see graphically (note though that since Propred uses 9mers and NetMHC 15mers, there's a shift)

../score.py --fa in/2sak_A.fasta --plot_hits_file out/2sak_A.pp5.pdf
../score.py --fa in/2sak_A.fasta --netmhcii --epi_thresh 10 --allele_set greenbaum11 --plot_hits_file out/2sak_A.nm10g.pdf
../score.py --fa in/2sak_A.fasta --netmhcii --epi_thresh 10 --allele_set paul15 --plot_hits_file out/2sak_A.nm10p.pdf

#7 Just the strong (2%) binders for the NetMHCII Greenbaum set

../score.py --fa in/2sak_A.fasta --netmhcii --epi_thresh 2 --allele_set greenbaum11 --plot_hits_file out/2sak_A.nm2g.pdf

#8 Individual scores for peptides, using predefined file of 9mers from SakSTAR

../score.py --pep in/2sak_A.pep > out/2sak_A_pep.txt

#9 Passing directly via command line

../score.py --report hits YLMVNVTGV

# 1,YLMVNVTGV,5,1,4,3,1,-,2,-,-

#10 Or from stdin

../score.py < in/2sak_A.seq

# * seq0 54

#11 Or from pdb file

../score.py --pdb in/2sak_A.pdb

# * 2sak_A_A 54

#12 Total scores for each of the three

../score.py --fsa in/all.fasta

# * 1eer_A 66
# * 2sak_A 54
# * 3r2x_C 37

#13 Hit reports and plots for each of the three proteins, using Propred 5% and NetMHCII Greenbaum 10%
# Use the $ wildcard to tell the script to substitute the appropriate sequence's name in the output filename

../score.py --fsa in/all.fasta --report hits --report_file out/$.pp.csv --plot_hits_file out/$.pp5.pdf

../score.py --fsa in/all.fasta --report hits --netmhcii --epi_thresh 10 --allele_set greenbaum11 --report_file out/$.nm10g.csv --plot_hits_file out/$.nm10g.pdf

#TODO: a file with a bunch of designed variants of a single protein to be scored for comparison

# ======================================== 
# Database creation

# There's really no reason to create a database from Propred, since it's very fast anyway and can be directly invoked from within Rosetta. However, since it is faster and its peptides are shorter, it may be useful to try it first on your protein, as NetMHCII will be slower and have more peptides due to the combinatoric difference between 15mer vs. 9mer. You can also use the positions specification and an aggressive pssm filter to keep things under control.

#1 Just wild-type epitopes, now stashed in a database

../db.py --fa in/2sak_A.fasta out/2sak_A.wt.pp5.db
sqlite3 out/2sak_A.wt.pp5.db 'select peptide,score from epitopes where score>0;' > out/2sak_A.hits.pp5.txt

#2 The database can then be used for scoring

../score.py --fa in/2sak_A.fasta --db out/2sak_A.wt.pp5.db

# * 2sak_A 54.0

#3 Set up a case where an epitope is missing from the database

cp out/2sak_A.wt.pp5.db out/2sak_A.wt_missing.pp5.db
sqlite3 out/2sak_A.wt_missing.pp5.db 'delete from epitopes where peptide="ALDATAYKE";'
../score.py --fa in/2sak_A.fasta --db out/2sak_A.wt_missing.pp5.db 

# * 2sak_A 154.0

../score.py --fa in/2sak_A.fasta --db out/2sak_A.wt_missing.pp5.db --db_unseen error

# ERROR unscored peptide ALDATAYKE

#4 A simple csv file allowing some mutations in the region in Choi Fig. 6 (native numbering). Save out the peptides for convenient inspection

../db.py --fa in/2sak_A_native.fasta --aa_csv in/2sak_A.region_muts.csv --peps_out out/2sak_A.region.pep out/2sak_A.region.pp5.db

# While the peptides aren't ordered (a set has been used to eliminate possible duplicates), you can see for example in out/2sak_A.region.pep the peptide at 67 ALDATAYKE along with substitutions of S for the first A and/or Q for the third A as allowed by the CSV file (=> SLDATAYKE, ALDATQYKE, SLDATQYKE), etc.

# Can also now ask for scores of sequences using these mutations

../score.py --db out/2sak_A.region.pp5.db ALDATAYKEFRVVELDPSAKI SLDATQYKEFRVVELDPSAKI ALDATAYKEFRTVELDPSAKI ALDATAYKEFRTTELDPSAKI

#5 Mutations allowed according to the pssm, but only in the region in Choi Fig. 6 and locked down as in Table 1. Positions are adjusted for 1-based pssm and Rosetta pose numbering (native-15 for SakSTAR and native-5 for HB36). For demo purposes (to keep it small and fast), raise the pssm threshold.

../db.py --fa in/2sak_A.fasta --pssm in/2sak_A.pssm --pssm_thresh 2 --positions 52-72 --lock 54,55,58,60,61 --res_out out/2sak_A.52-72.res --chain A out/2sak_A.52-72.pssm.pp5.db

../db.py --fa in/1eer_A.fasta --pssm in/1eer_A.pssm --pssm_thresh 2 --positions 140-160 --lock 150,151,155 --res_out out/1eer_A.140-160.res --chain A out/1eer_A.140-160.pssm.pp5.db

../db.py --fa in/3r2x_C.fasta --pssm in/3r2x_C.pssm --pssm_thresh 2 --positions 45-65 --lock 48,51,52,55,56,58,64 --res_out out/3r2x_C.45-65.res --chain C out/3r2x_C.45-65.pssm.pp5.db

#TODO: something interesting from this preprocessing just as-is, before heading into design
