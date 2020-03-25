#!/bin/bash

mkdir -p out

#1 List Propred epitope hits (5%) for SakStar, using the default Southwood set of 8 alleles

../mhc_score.py --fa in/2sak_A.fasta --report hits 

#2 Using native numbering, for ref against Choi Fig. 6

../mhc_score.py --fa in/2sak_A_native.fasta --report hits

#3 While the printout is useful from the command-line, the csv file is better for posterity and potentially for further analysis

../mhc_score.py --fa in/2sak_A.fasta --report hits --report_file out/2sak_A.pp5.csv

#4 Compare to NetMHCII (10% -- "weak binder" as of NetMHCII 2.3 web server), for the Greenbaum set of alleles

../mhc_score.py --fa in/2sak_A.fasta --report hits --netmhcii --epi_thresh 10 --allele_set greenbaum11 --report_file out/2sak_A.nm10g.csv

#5 And the Paul set

../mhc_score.py --fa in/2sak_A.fasta --report hits --netmhcii --epi_thresh 10 --allele_set paul15 --report_file out/2sak_A.nm10p.csv

#6 These things can be easier to see graphically (note though that since Propred uses 9mers and NetMHC 15mers, there's a shift)

../mhc_score.py --fa in/2sak_A.fasta --plot_hits_file out/2sak_A.pp5.pdf
../mhc_score.py --fa in/2sak_A.fasta --netmhcii --epi_thresh 10 --allele_set greenbaum11 --plot_hits_file out/2sak_A.nm10g.pdf
../mhc_score.py --fa in/2sak_A.fasta --netmhcii --epi_thresh 10 --allele_set paul15 --plot_hits_file out/2sak_A.nm10p.pdf

#7 Just the strong (2%) binders for the NetMHCII Greenbaum set

../mhc_score.py --fa in/2sak_A.fasta --netmhcii --epi_thresh 2 --allele_set greenbaum11 --plot_hits_file out/2sak_A.nm2g.pdf

#8 Individual scores for peptides, using predefined file of 9mers from SakSTAR

../mhc_score.py --pep in/2sak_A.pep > out/2sak_A_pep.txt

#9 Passing directly via command line

../mhc_score.py --report hits YLMVNVTGV

# 1,YLMVNVTGV,5,1,4,3,1,-,2,-,-

#10 Or from stdin

../mhc_score.py < in/2sak_A.seq

# * seq0 54

#11 Or from pdb file

../mhc_score.py --pdb in/2sak_A.pdb

# * 2sak_A_A 54

#12 Total scores for each of the three

../mhc_score.py --fsa in/all.fasta

# * 1eer_A 66
# * 2sak_A 54
# * 3r2x_C 37

#13 Hit reports and plots for each of the three proteins, using Propred 5% and NetMHCII Greenbaum 10%
# Use the $ wildcard to tell the script to substitute the appropriate sequence's name in the output filename

../mhc_score.py --fsa in/all.fasta --report hits --report_file out/$.batch.pp5.csv --plot_hits_file out/$.batch.pp5.pdf

../mhc_score.py --fsa in/all.fasta --report hits --netmhcii --epi_thresh 10 --allele_set greenbaum11 --report_file out/$.batch.nm10g.csv --plot_hits_file out/$.batch.nm10g.pdf

#TODO: a file with a bunch of designed variants of a single protein to be scored for comparison

