This file contains one example each of the two supported PSSM formats for use with mhc_gen_db.py.

mhc_gen_db.py has been implemented with some flexibility to allow for slightly different PSSM formats.  The results 
should be carefully checked if you try an unsupported format.

In general, mhc_gen_db.py will search for a AA order line, which contains the 20 AAs (whitespace delimited).  The command 
line psiblast output actually contains two copies of this, as it contains two tables, one with the log-transformed 
probabilities of encountering a specific AA, and one with percentages.  The ASN-to-matrix PSSM converter, for use 
with the online psiblast outputs the AA header row starting with P (position), then C (consensus), and finally the 
20 AAs.

It is critical that after the AA header line, the whitespace-delimited columns consist of (1) a protein position 
number (corresponding to the input FASTA/PDB file), (2) the "consensus/wild-type" AA identity, and (3-22) the 
scores for the 20 AAs, given as "log-odds" scores, which are the log(base 2) of the observed substitution frequency 
at a given position divided by the expected substitution frequency at that position.

command_line_pssm.txt is the pssm produced by the command line psiblast tool.
online_pssm.txt is the pssm produced using the online psiblast tool and converted by https://www.ncbi.nlm.nih.gov/Class/Structure/pssm/pssm_viewer.cgi