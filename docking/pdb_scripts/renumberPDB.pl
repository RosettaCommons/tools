#!/usr/bin/perl

# renumberPDB.pl: renumber the residue numbers in a PDB file
#
# this script is for cleaning files by removing insert letters,
# particularly for loading antibody pdbs into Tanja's code
#
# numbering is started with the first residue in a chain, and
# incremented by ones from there
#
# beware, one will lose track of missing residues, as discontinuous 
# numbering now becomes continuous



# JJG 12/11/01



if ($#ARGV < 0) {
	print STDERR "usage: $0 <pdbfile>\n";
	print STDERR "sequentially renumber the residues in the pdb file\n";
	exit -1;
}

$pdbfile = shift @ARGV;
open (PDB, $pdbfile);

$n=-1; # current residue number, indicate new chain with -1
while (<PDB>) {

    $line = $_;
    $last_chain = $chain;
    $chain = substr($line,21,1);
    if ($chain ne $last_chain) {$n=-1;}  # start over if we find a new chain
    if (/TER/) {$n=-1;}  # start over if we find a terminus
	
    if (/ATOM/) {
	if ($n == -1) { # start with the starting number in the chain
	    $n = substr($line,22,4);
	    $last_res_id = substr($line,22,5);
	}
	$res_id = substr($line,22,5);
	if ($res_id ne $last_res_id) {
	    $n++;
	    $last_res_id = $res_id;
	}
	substr($line,22,5) = sprintf "%4d ",$n;  
    }
    print $line;
}
