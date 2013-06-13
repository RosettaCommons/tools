#!/usr/bin/perl

# pdb_subtract_scores.pl : compare two pdbs by subracting the scores 
# contained within them (ie scores from rosetta)



if ($#ARGV < 1) {
	print STDERR "usage: $0 <pdbfile1> <pdbfile2> \n";
	print STDERR "Compare scores residue by residue\n";
	exit -1;
}
$pdbf1 = shift @ARGV;
$pdbf2 = shift @ARGV;

open (PDB1, $pdbf1);
open (PDB2, $pdbf2);

#skip the coordinates
while (<PDB1>) {  last if ( /^TER/ ); }
while (<PDB1>) {  last if ( /^TER/ ); }

while (<PDB2>) {  last if ( /^TER/ ); }
while (<PDB2>) {  last if ( /^TER/ ); }


while (!eof(PDB1)) {

    $line1 = <PDB1>;
    $line2 = <PDB2>;

    @fields1 = split(" ",$line1);
    @fields2 = split(" ",$line2);

#    print "F1: @fields1";
#    print "F2: @fields2";

    print "$fields1[0]\t";
    if ( $fields1[0] != $fields2[0] ) {
	printf "\nMismatch! \n $line1 \n $line2 \n";
	exit;
    };
    
    for ($i=1; $i<scalar(@fields1); $i++) {
	if ( $fields1[$i] =~ /\d/) {
	    #printf "Y$fields2[$i]-$fields1[$i]=";
	    $newf[$i] = $fields2[$i] - $fields1[$i];
	    $newf[$i] = int(100*$newf[$i])/100;
	}
	else {
	    $newf[$i] = $fields2[$i];
	}
	printf "$newf[$i]\t";
    }
    printf "\n";
}




