The examples here are based on the case studies used in an initial Rosetta deimmunization method by Y. Choi, K. Griswold, and C. Bailey-Kellogg, "Structure-based Redesign of Proteins for Minimal T-cell Epitope Content",  J. Comp. Chem 2013, pmid 23299435. The mhc_epitope_energy implementation is based on code from that publication, but brought up to date and extended and generalized. Thus the results will not match exactly.

The pdb files have been preprocessed by ROSETTA/tools/protein_tools/scripts/clean_pdb.py and thus all have "pose numbering" starting at 1. The corresponding fasta files likewise are assumed to start at 1, so that preprocessing for Rosetta will correspond. Note that the native start for SakSTAR (2sak) is 16 (so add 15 to pose numbers) and for HB-36 (3r2x_C) is 6 (so add 5). This is also important when "locking down" positions in the epitope database examples, not allowing substitutions at positions listed as fixed in Choi et al. table 1. 

For human-interpretable preprocessing in cases where the pdb file doesn't actually start at residue 1, the scripts allow for fasta files to specify a starting position and for pdb files to start wherever. The "native.fasta" files use that extension, in order to allow direct comparison to figure 6 in the Choi et al. paper (and results from immunogenicity studies). 

The pssm files were produced for the 2013 paper, so of course results from searches on the current database will differ. Command-line PSI-BLAST was used to generate a pssm directly. The web version generates an asn file that can be converted by a web tool (as of this writing, at https://www.ncbi.nlm.nih.gov/Class/Structure/pssm/pssm_viewer.cgi). 

The demo-score.sh file and demo-db.sh files have some examples of how to use the scripts, using inputs in "in/" and saving outputs in "out/"; reference outputs are saved in "ref-out/". For the sake of demonstration, in addition to the real fasta files,"2sak_A.pep" lists all overlapping 9mers from SakSTAR and "2sak_A.seq" has the sequence without the header, while "all.fsa" includes all three proteins.

The tools for dealing with experimental (from IEDB) are not as far along, but demo-expt.sh illustrates getting data and formatting into files that can be used in scoring and design, the same as an epitope predictor.

The script demo-diff.sh will compare the relevant files from the out directory to the reference ref-out directory.  You can think of this as a bit like an integration test.
