Contact Steven Lewis, smlewi@gmail.com, for support.  Will Hansen wrote much of the script used in the last step.

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

"pdb*ent.gz" is a pattern that matches how the files are stored, it's not xxxx.pdb.gz as you'd expect.  Find will return 116000+ lines like:

/home/smlewis/Desktop/PDB_mirror/whole_PDB_as_PDBgz/is/pdb1iso.ent.gz
/home/smlewis/Desktop/PDB_mirror/whole_PDB_as_PDBgz/is/pdb4is9.ent.gz
/home/smlewis/Desktop/PDB_mirror/whole_PDB_as_PDBgz/is/pdb4ist.ent.gz
/home/smlewis/Desktop/PDB_mirror/whole_PDB_as_PDBgz/is/pdb3is8.ent.gz
/home/smlewis/Desktop/PDB_mirror/whole_PDB_as_PDBgz/is/pdb4isl.ent.gz
/home/smlewis/Desktop/PDB_mirror/whole_PDB_as_PDBgz/is/pdb2isw.ent.gz

Now allpdbs.l is a list file of all your inputs.

Next we randomize it.  This is optional, but I think it's a good idea, because I think there are stretches with a bunch of huge ribosome structures that are adjacent in PDB-ID space.  If you are going to run on more than one processor, you'll want them distributed; random is close enough.

sort -R  allpdbs.l > allpdbs.random.l

Next we turn our PDB target paths into Rosetta commands.  You can do this with rectangle edit in emacs (c-space for mark, c-x r t for rectangle insert), or sed.

From
/home/smlewis/whole_PDB_as_PDBgz/is/pdb1iso.ent.gz
/home/smlewis/whole_PDB_as_PDBgz/is/pdb4is9.ent.gz
/home/smlewis/whole_PDB_as_PDBgz/is/pdb4ist.ent.gz

We want
../../PDB_diagnostic.linuxgccrelease @../../options -s /home/smlewis/whole_PDB_as_PDBgz/el/pdb5els.ent.gz >>log 2>&1
../../PDB_diagnostic.linuxgccrelease @../../options -s /home/smlewis/whole_PDB_as_PDBgz/fa/pdb4fa8.ent.gz >>log 2>&1
../../PDB_diagnostic.linuxgccrelease @../../options -s /home/smlewis/whole_PDB_as_PDBgz/zk/pdb1zk1.ent.gz >>log 2>&1

(yes, the PDB codes changed, ignore that).  We've added a Rosetta invocation and options to the front end, and logging to the back end (this is for bash, BTW).

Sed command is:

cat allpdbs.random.l | sed 's:^:../../PDB_diagnostic.linuxgccrelease @../../options -s :' | sed 's:$: >>log 2>\&1:' > rosetta_commands.l

Note the escaping of one of the &s in sed.

#JD0
Go read and/or perform the JD0 readme.

#options
Let me suggest you use the options file in the PDB_diagnostic integration test as inspiration.  Don't use it verbatim (you'll want at least the pdb ccd library).

#Process the results
I tend to get on the order of 20 GB of log files when running this experiment, but the vast majority of those logs are certain output lines repeated many, many times inside (code) loops caused by (chemical) ring structures that the reading-in machinery has difficulty with.  You'll want to grep out some of these uninformative lines from the raw logs to speed later processing and slash the memory cap needed for the last step.  At a minimum, I suggest (note subdirs comes from JD0)

cat subdirs/*/log | grep -v "missing heavyatom" | grep -v " atoms at position " > slimlog

or (note -h to grep when using find)

find subdirs -name "log" | xargs grep -vh "missing....

Continue cutting slimlog as necessary to make the file size reasonable by removing repeated uninformative log lines; use your best judgment.

#Analyze the results

check_logs.py slimlog > check_logs.log

This will also spawn dozens of extra files.  The script parses the log file to determine which PDBs had errors, and can classify more than 20 of those errors.

check_logs.log is a summary of results.

*.list is just PDBs that failed with a particular error type.
*.cmdpath is the command paths (as in rosetta_commands.l above) ready to be resubmitted via JD0.
*.log is log snippets for explanation and diagnosis.
