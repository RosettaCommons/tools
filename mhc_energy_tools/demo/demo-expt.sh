mkdir -p out

# Download a fresh snapshot of IEDB and load it into a local database
# For this, need to have mysql / mariadb up and running; example assumes login as root / no password
# Also need to have installed mysql.connector https://dev.mysql.com/downloads/connector/python/
# The download can take a while, so this example stops after the download and then the next example does something with it; it's possible also to do both in a single execution

# I've left it commented out so you have to really want to execute it to do so, by copy-paste
#../mhc_data_db.py --iedb_fresh_mysql iedb

# Pull out the DRB1*01:01 data (the 'test' allele-set), store all 9mer cores as potential epitopes
../mhc_data_db.py --iedb_mysql iedb --allele_set test --csv out/iedb_0101_allcores.csv

# Or just the good-enough cores according to Propred
# Relies on 'test' allele-set being the same for both
../mhc_data_db.py --iedb_mysql iedb --allele_set test --propred --cores predicted_good --csv out/iedb_0101_pp5cores.csv

# Now can scan a protein with this as the scorer
# Note that the unseen score here is set to 0, since not being seen is good (the default is a penalty)
../mhc_score.py --fa in/2sak_A.fasta --csv out/iedb_0101_allcores.csv --db_unseen score --db_unseen_score 0

# 2sak_A 0

# Check this one out! Epo has been completely mapped, and apparently every peptide is a hit.
../mhc_score.py --fa in/1eer_A.fasta --csv out/iedb_0101_allcores.csv --db_unseen score --db_unseen_score 0

# 1eer_A 109.0

# Do the same with a file saved from IEDB. Since 1eer turned out to be interesting, I downloaded its MHC II ligand assay data (from assay tab)

../mhc_data_db.py --iedb_csv in/iedb_epo_mhcii_assays.csv --allele_set test --csv out/iedb_epo_0101_allcores.csv
../mhc_score.py --fa in/1eer_A.fasta --csv out/iedb_epo_0101_allcores.csv --db_unseen score --db_unseen_score 0

# 1eer_A 109.0
