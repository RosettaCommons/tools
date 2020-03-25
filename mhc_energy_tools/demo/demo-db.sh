#!/bin/bash

mkdir -p out

# There's really no reason to create a database from Propred, since it's very fast anyway and can be directly invoked from within Rosetta. However, since it is faster and its peptides are shorter, it may be useful to try it first on your protein, as NetMHCII will be slower and have more peptides due to the combinatoric difference between 15mer vs. 9mer. You can also use the positions specification and an aggressive pssm filter to keep things under control.

#1 Just wild-type epitopes, now stashed in a database

../mhc_gen_db.py --propred --fa in/2sak_A.fasta --db out/2sak_A.wt.pp5.db
sqlite3 out/2sak_A.wt.pp5.db 'select peptide,score from epitopes where score>0;' > out/2sak_A.hits.pp5.txt

#1b Can also store directly to csv

../mhc_gen_db.py --propred --fa in/2sak_A.fasta --csv out/2sak_A.wt.pp5.csv
# Sort the csv file to make the diff pass, leaving the header line alone
(read -r; printf "%s\n" "$REPLY"; sort) < out/2sak_A.wt.pp5.csv > out/tmp
mv out/tmp out/2sak_A.wt.pp5.csv

#2 The database can then be used for scoring

../mhc_score.py --fa in/2sak_A.fasta --db out/2sak_A.wt.pp5.db

# * 2sak_A 54.0

#2b And so can the csv

../mhc_score.py --fa in/2sak_A.fasta --csv out/2sak_A.wt.pp5.csv

# * 2sak_A 54.0

#3 Can convert back and forth

../mhc_gen_db.py --csv_in out/2sak_A.wt.pp5.csv --db out/2sak_A.wt.pp5.v2.db
../mhc_score.py --fa in/2sak_A.fasta --db out/2sak_A.wt.pp5.v2.db

# * 2sak_A 54.0

../mhc_gen_db.py --db_in out/2sak_A.wt.pp5.v2.db --csv out/2sak_A.wt.pp5.v2.csv
# Sort the csv file to make the diff pass, leaving the header line alone
(read -r; printf "%s\n" "$REPLY"; sort) < out/2sak_A.wt.pp5.v2.csv > out/tmp
mv out/tmp out/2sak_A.wt.pp5.v2.csv
../mhc_score.py --fa in/2sak_A.fasta --csv out/2sak_A.wt.pp5.v2.csv

# * 2sak_A 54.0

#4 Set up a case where an epitope is missing from the database

cp out/2sak_A.wt.pp5.db out/2sak_A.wt_missing.pp5.db
sqlite3 out/2sak_A.wt_missing.pp5.db 'delete from epitopes where peptide="ALDATAYKE";'
../mhc_score.py --fa in/2sak_A.fasta --db out/2sak_A.wt_missing.pp5.db 

# * 2sak_A 154.0

../mhc_score.py --fa in/2sak_A.fasta --db out/2sak_A.wt_missing.pp5.db --db_unseen error

# ERROR unscored peptide ALDATAYKE

#5 A simple csv file allowing some mutations in the region in Choi Fig. 6 (native numbering). Save out the peptides for convenient inspection

../mhc_gen_db.py --propred --fa in/2sak_A_native.fasta --aa_csv in/2sak_A.region_muts.csv --peps_out out/2sak_A.region.pep --db out/2sak_A.region.pp5.db

# While the peptides aren't ordered (a set has been used to eliminate possible duplicates), you can see for example in out/2sak_A.region.pep the peptide at 67 ALDATAYKE along with substitutions of S for the first A and/or Q for the third A as allowed by the CSV file (=> SLDATAYKE, ALDATQYKE, SLDATQYKE), etc.

# Can also now ask for scores of sequences using these mutations

../mhc_score.py --db out/2sak_A.region.pp5.db ALDATAYKEFRVVELDPSAKI SLDATQYKEFRVVELDPSAKI ALDATAYKEFRTVELDPSAKI ALDATAYKEFRTTELDPSAKI

#6 Mutations allowed according to the pssm, but only in the region in Choi Fig. 6 and locked down as in Table 1. Positions are adjusted for 1-based pssm and Rosetta pose numbering (native-15 for SakSTAR and native-5 for HB36). For demo purposes (to keep it small and fast), raise the pssm threshold.

../mhc_gen_db.py --propred --fa in/2sak_A.fasta --pssm in/2sak_A.pssm --pssm_thresh 2 --positions 52-72 --lock 54,55,58,60,61 --res_out out/2sak_A.52-72.res --chain A --db out/2sak_A.52-72.pssm.pp5.db

../mhc_gen_db.py --propred --fa in/1eer_A.fasta --pssm in/1eer_A.pssm --pssm_thresh 2 --positions 140-160 --lock 150,151,155 --res_out out/1eer_A.140-160.res --chain A --db out/1eer_A.140-160.pssm.pp5.db

../mhc_gen_db.py --propred --fa in/3r2x_C.fasta --pssm in/3r2x_C.pssm --pssm_thresh 2 --positions 45-65 --lock 48,51,52,55,56,58,64 --res_out out/3r2x_C.45-65.res --chain C --db out/3r2x_C.45-65.pssm.pp5.db

#TODO: something interesting from this preprocessing just as-is, before heading into design
