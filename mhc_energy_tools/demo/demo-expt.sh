#!/bin/bash

# This demo runs a series of tests that use the IEDB download capabilities.
# For this, need to have mysql / mariadb up and running, with access to a database.
# Also need to have installed mysql.connector https://dev.mysql.com/downloads/connector/python/
# By default, we assume that the database is called 'iedb' and you have root access with no password.
# Arguments to the script will allow you to specify a database name, username, and password.

#argparser for mysql settings
dbname=iedb #Default mysql database name
while getopts ":d:u:p:" opt; do
	case ${opt} in
		d ) dbname=$OPTARG
			;;
		u ) uname=$OPTARG
			;;
		p ) pword=$OPTARG
			;;
		\? ) echo "Usage: demo-expt.sh [-d mysql_database_name] [-u mysql_database_userid] [-p mysql_password]"
			echo "The default mysql_database_name is iedb.  By default, no username or password are provided."
			exit 1
			;;
		: ) echo "Invalid option: $OPTARG requires an argument, if provided."
			exit 1
			;;
	esac
done

#Make the "sql string"
#Ensures that we only pass the options with values to the python script
if [ $uname ]; then
	sqlstring="$sqlstring --mysql_user $uname"
fi
if [ $pword ]; then
	sqlstring="$sqlstring --mysql_pw $pword"
fi

mkdir -p out

# Download a fresh snapshot of IEDB and load it into a local database
# I've left it commented out so you have to really want to execute it to do so, by copy-paste
# Subsequent commands will not work unless you already have downloaded the db using this command, or you uncomment it.
#../mhc_data_db.py --iedb_fresh_mysql iedb
# With arguments:
#../mhc_data_db.py --iedb_fresh_mysql iedb --mysql_user sqlusername --mysql_pw sqlpassword
# Uncomment and run in this script:
#../mhc_data_db.py --iedb_fresh_mysql $dbname $sqlstring

# Pull out the DRB1*01:01 data (the 'test' allele-set), store all 9mer cores as potential epitopes
../mhc_data_db.py --iedb_mysql $dbname $sqlstring --allele_set test --csv out/iedb_0101_allcores.csv
# Sort the csv file to make the diff pass, leaving the header line alone
(read -r; printf "%s\n" "$REPLY"; sort) < out/iedb_0101_allcores.csv > out/tmp
mv out/tmp out/iedb_0101_allcores.csv

# Or just the good-enough cores according to Propred
# Relies on 'southwood98' allele-set being the same for both
../mhc_data_db.py --iedb_mysql $dbname $sqlstring --allele_set southwood98 --propred --cores predicted_good --csv out/iedb_0101_pp5cores.csv
# Sort the csv file to make the diff pass, leaving the header line alone
(read -r; printf "%s\n" "$REPLY"; sort) < out/iedb_0101_pp5cores.csv > out/tmp
mv out/tmp out/iedb_0101_pp5cores.csv

# Now can scan a protein with this as the scorer
# Note that the unseen score here is set to 0, since not being seen is good (the default is a penalty)
../mhc_score.py --fa in/2sak_A.fasta --csv out/iedb_0101_allcores.csv --db_unseen score --db_unseen_score 0

# 2sak_A 0

# Check this one out! Epo has been completely mapped, and apparently every peptide is a hit.
../mhc_score.py --fa in/1eer_A.fasta --csv out/iedb_0101_allcores.csv --db_unseen score --db_unseen_score 0

# 1eer_A 109.0

# Do the same with a file saved from IEDB. Since 1eer turned out to be interesting, I downloaded its MHC II ligand assay data (from assay tab)

../mhc_data_db.py --iedb_csv in/iedb_epo_mhcii_assays.csv --allele_set test --csv out/iedb_epo_0101_allcores.csv
# Sort the csv file to make the diff pass, leaving the header line alone
(read -r; printf "%s\n" "$REPLY"; sort) < out/iedb_epo_0101_allcores.csv > out/tmp
mv out/tmp out/iedb_epo_0101_allcores.csv
../mhc_score.py --fa in/1eer_A.fasta --csv out/iedb_epo_0101_allcores.csv --db_unseen score --db_unseen_score 0

# 1eer_A 109.0
