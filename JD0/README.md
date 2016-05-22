This is the hacky shell script job distributor Steven Lewis wrote at ChemXRW 2016 for the "try reading the whole PDB" experiments.
It's not robust in any sense of the term except that it lets you run many, many runs of Rosetta that can't crash each other the way a large -l run will.

I suggest also reading the chem_xrw README for an explanation of when you'd use this and some framing.  It will also explain why these instructions are steps 3 and 4.

#What you start with
On your own (or with chem_xrw), assemble a file with a ton of Rosetta command lines, like this, named rosetta_commands.list or similar

../../PDB_diagnostic.linuxgccrelease @../../options -s /home/smlewis/whole_PDB_as_PDB/el/pdb5els.ent.gz >>log 2>&1
../../PDB_diagnostic.linuxgccrelease @../../options -s /home/smlewis/whole_PDB_as_PDB/fa/pdb4fa8.ent.gz >>log 2>&1
../../PDB_diagnostic.linuxgccrelease @../../options -s /home/smlewis/whole_PDB_as_PDB/zk/pdb1zk1.ent.gz >>log 2>&1
../../PDB_diagnostic.linuxgccrelease @../../options -s /home/smlewis/whole_PDB_as_PDB/jb/pdb4jbd.ent.gz >>log 2>&1
../../PDB_diagnostic.linuxgccrelease @../../options -s /home/smlewis/whole_PDB_as_PDB/gy/pdb2gyx.ent.gz >>log 2>&1

#Step 3: setup

setup.bash has that file name hardcoded into its "split" command line (line 8).  Fix that hardcoding as necessary.  This script is written to make 10 subdirectories (to run on 10 processors).  You'll need to read up on split to learn how to recraft that command line for other numbers, but here are some guesses:

1-9, where # is the number you need
split -a 1 -n l/#  --numeric-suffixes=1 rosetta_commands.list allpdbs_slice_

10-99, same as it's saved in git, just change the 10

100-999, change the -a 2 to -a 3 and the 10 to 100 or 355 or whatever.

10 is also hardcoded into the for loop seq subcommand, change that as needed.

The script will make 10 subdirectories in a subdirs subdirectory, each with its own slice of the rosetta_commands.list, ready to be run. Notice that you embedded relative paths ../../ in rosetta_commands.txt: thus you just need one copy of your options file and Rosetta exe symlink, in the directory above subdirs.

#Step 4: run
submitter.bash must also be edited for the right number of processors (the same seq command as setup.bash).  It uses "nohup" here to just run the command without interruptions (still OK to kill or pause), but you'd substitute nohup for your cluster command submission as needed.