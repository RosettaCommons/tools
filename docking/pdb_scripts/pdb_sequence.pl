#!/usr/bin/perl

if ($#ARGV < 0) {
        print STDERR "usage: $0 <pdbfile> [simple] \n";
	print STDERR "print the single-letter sequence from a pdb\n";
	print STDERR "gaps in residue numbering indicated by -\n";
	print STDERR "and insertions (ie. residue 15A) indicated by +\n";
	print STDERR "simple option excludes insertion/missing residue indicators\n";
        exit -1;
}

$pdbfile=$ARGV[0];
$simple="true" if ($ARGV[1] eq "simple");

@residues=`grep ^ATOM $pdbfile | cut -c 18-27 |uniq`;
#print @residues;
chomp(@residues);

$nres=0;
$i=0;
for $r (@residues) {
    $three = substr($r,0,3);
    $ch = substr($r,4,1);
    $n = substr($r,5,4);
#    print "ch $ch n $n --$r";

    if (@seq3 && @seq1 && $ch ne $last_ch) {
	print "---$last_ch---\n";
	#print "@seq3\n";
	while (@seq1){ print splice(@seq1,0,60),"\n"; }
	undef @seq3;
	$nres=0;
    }
    elsif (@seq3 && ($n != $last_n + 1) && !$simple) {
	if ( $n == $last_n ) {
	    push(@seq3,'+++');
	    push(@seq1,'+');
	}
	else {
	    push(@seq3,'---');
	    push(@seq1,'-');
	}
    }

    push(@seq3,$three);
    $one = &one_from_three($three);
    push(@seq1,$one);
    $last_ch = $ch;
    $last_n  = $n;
    $nres++;
}

print "---$last_ch---\n";
#print "@seq3\n";
while (@seq1){ print splice(@seq1,0,60),"\n"; }


exit 0;




sub one_from_three {
    local $three = shift;
    local $one = ".";
    $one="A" if ($three eq "ALA");
    $one="C" if ($three eq "CYS");
    $one="D" if ($three eq "ASP");
    $one="E" if ($three eq "GLU");
    $one="F" if ($three eq "PHE");
    $one="G" if ($three eq "GLY");
    $one="H" if ($three eq "HIS");
    $one="I" if ($three eq "ILE");
    $one="K" if ($three eq "LYS");
    $one="L" if ($three eq "LEU");
    $one="M" if ($three eq "MET");
    $one="N" if ($three eq "ASN");
    $one="P" if ($three eq "PRO");
    $one="Q" if ($three eq "GLN");
    $one="R" if ($three eq "ARG");
    $one="S" if ($three eq "SER");
    $one="T" if ($three eq "THR");
    $one="V" if ($three eq "VAL");
    $one="W" if ($three eq "TRP");
    $one="Y" if ($three eq "TYR");
    return $one;
}
