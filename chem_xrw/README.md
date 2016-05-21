This directory describes the process of running the PDB_diagnostic app against the whole PDB.  PDB_diagnostic was created for the purpose of exercising a the code to generate Poses from PDB-format structure files (pdbs and mmcifs), and is a pilot app over in the main repo.

The steps are:
1) Get the whole PDB (some variant of rsyncPDB.sh)
2) Generate Rosetta command lines, one for each input structure
3-4) Generate folders in which to run these jobs, and run them (this is covered in the companion JD0 folder)
5) Process the results
6) Analyze the results (check_logs.py in this folder).

# Get the whole PDB
The PDB maintains an excellent rsync script described at http://www.rcsb.org/pdb/static.do?p=general_information/news_publications/newsletters/2003q3/focus_rsync.html and hosted at ftp://ftp.rcsb.org/pub/pdb/software/ .  It will explain how to download a mirror of the PDB.  Do note the odd structure, where files are in folders that match the MIDDLE TWO characters of their PDB ids.  This tends to keep codeposited (contiguously-numbered) structures together, while balancing the number of folders vs. files per folder (inode limits, you know).

Note that as of this writing, mmCIF input is gzipped in the rsync mirror, but the library we use to read mmCIF expects gunzipped files, so you'll need to gunzip them all.  I suggest making a deep copy and then using:

find /path/to/copy -name "*cif.gz" -exec gunzip {} \;

Note also that this is NOT necessary for PDBs.

# Generate Rosetta command lines for each input
PDB_diagnostic is used to examine whether structures are read in successfully, and if not, why not.  A lot crash, and some even segfault.  You can't just run Rosetta in MPI with -l the_whole_pdb.list, the segfaults will kill it.  Thus, this tool uses JD0 to route around the problem.  JD0 means "no job distribution at all" - you have to do it manually with one Rosetta invocation per input.

Here's how I do it.  First, generate a file containing the paths to every file in your input.  Find is the easy way:

find /path/to/pdbmirror -name "pdb*ent.gz" > allpdbs.l

"pdb*ent.gz" is a pattern t






This script check_logs.py, along with the Rosetta application PDB_diagnostic, is used to determine what PDBs Rosetta has problems reading.  PDB_diagnostic is intended to be run against the entire PDB dataset (all 116K+ PDBs).  The log files are then concatenated and fed through check_logs.py.  The script parses the log file to determine which PDBs had errors, and can classify more than 20 of those errors.  Contact Steven Lewis, smlewi@gmail.com, for support.  Will Hansen wrote much of the script.
