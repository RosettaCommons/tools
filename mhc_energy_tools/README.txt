Rosetta mhc_epitope_energy tools, developed by Chris Bailey-Kellogg (cbk@cs.dartmouth.edu) and Brahm Yachnin (brahm.yachnin@rutgers.edu)
Initial release, Fall 2018; updated June 2019

These tools support the use of the Rosetta mhc_epitope_energy term in redesigning proteins so as to reduce their epitope content. See documentation on that energy term for further background on using it in a design process. The energy term incorporates scores for the constituent peptides of a protein, by either invoking a matrix-based epitope predictor or querying a database. Matrix-based predictors are fast and may be directly invoked during the optimization. A database allows ready incorporation of known epitope data. A database also enables incorporation of external predictors, such as NetMHCII, that only provide a command-line interface and thus are not amenable to using during design. This type of database can be generated using mhc_gen_db.py. The general approach is to up front (before performing design) precompute and store scores of possible peptides that may be considered. This requires limiting the space of possible mutations, but that can be beneficial to focus the design anyway. Alternatively, a database based on experimentally validated epitopes can be generated using mhc_data_db.py.

In general, a "database" can be either a bonafide sqlite3 database or just a csv file. Since python and Rosetta both support sqlite3, there's no barrier to using that -- it's fairly lightweight and random access. However, a csv file is easier to look at manually and perhaps to manipulate with other programs.

The python tools help to preprocess a wild-type sequence and possible mutations before performing design, as well as to postprocess possible designs.

* Scoring: mhc_score.py. Provide sequence(s) and specify a predictor and its parameters; the script generates either the total score or a peptide-by-peptide report or plot. Numerous command-line arguments allow control over the options; invoke "score.py --help" for a description.

* Database creation: mhc_gen_db.py. Provide a sequence and specification of possible mutations; the script generates and scores epitopes and creates/augments the database with them. 

* Experimental database creation: mhc_data_db.py. Provide experimental information (currently IEDB format csv or mysql database); the script organizes into a database for use in scoring and design.

Note that mhc_score.py can either directly invoke a predictor or use a database.

A few notes about installation:
* These were developed for python 3
* To find the matrices for a matrix-based predictor (e.g., Propred), the script looks where it thinks the Rosetta database should be, assuming that the tools are installed in the same relative directory structure. It also checks in $ROSETTA/main/database/scoring/score_functions/mhc_epitope, if environment variable $ROSETTA is set.
* A simple interface is provided to NetMHCII, version 2.3. That software must be installed separately, and the environment variable $NMHOME set to point to the directory holding the netMHC2-3 script.

More details about these scripts can be found in documentation/rosetta_basics/scoring/mhc-energy-tools (and on the Rosetta Wiki).
