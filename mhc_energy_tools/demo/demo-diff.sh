#!/bin/bash
set -e #If a command fails, stop

#This script diffs the output from the demos.
#log files (especially demo-expt.log and demo-score.log) cannot be compared because of tmp directories and random order of reading some of the db data.

echo Preparing to compare all demo output to the \"ref-out\" directory.
echo Note that if you haven\'t run all of the demos, some of the files will be listed as \"missing\".
echo Also note that the output from demo-expt.sh may differ as the IEDB.org database is updated.  For that matter, updates to NetMHCII may also change things.
echo Press Enter to continue.
read
echo

echo Comparing text files \(.csv, .txt, .res, .pep\)...
#diff the two directories, excluding pdf, log and db files.
#To prevent silent output, uses -s to name identical files as well as differing files.
diff out ref-out -s -x '*.pdf' -x '*.log' -x '*.db'
echo Done comparing text files.  Press Enter to continue.
read
echo

echo Comparing SQL .db files...
mkdir tmpsql
mkdir ref-tmpsql

#Loop over .db files in out
for db in out/*.db; do
	#The filename without an extension
	base=`basename $db`
	base=${base%%.db}
	
	#Each db file should have alleles, epitopes, and meta tables
	#For each db file, dump contents of all three tables into a text file in tmpsql (for out .db files) or ref-tmpsql (for ref-out)
	#The data need to be sorted, as the order of the peptides in the table isn't enforced.
	sqlite3 $db '.dump alleles' | sort > tmpsql/`basename $db`
	sqlite3 ref-out/`basename $db` '.dump alleles' | sort > ref-tmpsql/`basename $db`
	sqlite3 $db '.dump epitopes' | sort >> tmpsql/`basename $db`
	sqlite3 ref-out/`basename $db` '.dump epitopes' | sort >> ref-tmpsql/`basename $db`
	sqlite3 $db '.dump meta' | sort >> tmpsql/`basename $db`
	sqlite3 ref-out/`basename $db` '.dump meta' | sort >> ref-tmpsql/`basename $db`
done

#Now compare the text files in tmpsql and ref-tmpsql, which conveniently have the same names as the original database files
diff tmpsql ref-tmpsql -s
echo Done comparing SQL .db files.
echo The pdf files and log files cannot be compared automatically, but if everything else is OK, those should be too.

#Cleanup carefully
rm tmpsql/*.db
rm ref-tmpsql/*.db
rmdir tmpsql
rmdir ref-tmpsql
