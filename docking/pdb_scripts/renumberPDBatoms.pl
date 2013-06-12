#!/usr/bin/perl

# renumberPDB.pl: renumber the atom numbers in a PDB file
#
# this script is for making atom # continuous, so as to avoid
#confusing rasmol and other programs that require continuous numbering.
#
# numbering is started with the first atom in a chain, and
# incremented by ones from there
#
# beware, one will lose track of missing atoms, as discontinuous 
# numbering now becomes continuous
#
#This script is a slightly modified version of renumberPDB.pl


# JJG 12/11/01
#Mike Daily 2/26/03

#a test comment for CVS

if ($#ARGV < 0) {
	print STDERR "usage: $0 <pdbfile>\n";
	print STDERR "sequentially renumber the atoms in the pdb file\n";
	exit -1;
}

$pdbfile = shift @ARGV;
open (PDB, $pdbfile);

$n=1; # current atom number
while (<PDB>) {

    $line = $_;
    if (/ATOM|HETATM/) {
	substr($line,7,5) = sprintf "%4d ",$n;  
	$n++;
    }
    print $line;
}
