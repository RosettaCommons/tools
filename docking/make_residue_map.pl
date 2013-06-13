#!/usr/bin/perl

if (@ARGV<3) {
    print STDERR "USAGE: $0 pdb1 pdb2 clustal.aln\n\n";
    print STDERR "Use the clustalw alignment to create a mapping of residue\n";
    print STDERR "identifiers between the two pdb files\n";
    exit;
}

$pdbfile1 = shift @ARGV;
$pdbfile2 = shift @ARGV;
$alnfile  = shift @ARGV;

#open(PDB1,$pdbfile1);
#open(PDB2,$pdbfile2);


#----read alignment
open(ALN,$alnfile);
$buf = <ALN>; # header line
$buf = <ALN>; # blank

while (!eof(ALN)) {
    $buf = <ALN>; # blank
    ($buf1 = <ALN>) =~ /\w+\s+(\S+)/; $buf1 = $1;
    ($buf2 = <ALN>) =~ /\w+\s+(\S+)/; $buf2 = $1;
    $buf = <ALN>; # blank

#    print "1: $buf1\n";
#    print "2: $buf2\n";

    push(@seq1,unpack("A1" x length($buf1), $buf1));
    push(@seq2,unpack("A1" x length($buf2), $buf2));
}

#print "1: @seq1\n";
#print "2: @seq2\n";


#----get residue lists from pdbs
@tmp = `grep ^ATOM $pdbfile1 | cut -c 17-27 |uniq`;
foreach (@tmp) {
    push(@reslist1,substr($_,5,6) . substr($_,0,4));
}
@tmp = `grep ^ATOM $pdbfile2 | cut -c 17-27 |uniq`;
foreach (@tmp) {
    push(@reslist2,substr($_,5,6) . substr($_,0,4));
}


#----map 'em
for ($i=0,$j1=0,$j2=0; $i<scalar(@seq1); $i++) {
    #chu: don't include gaps in both sequences.
    next if ($seq1[$i] eq '-' && $seq2[$i] eq '-');

    if ($seq1[$i] eq '-') { # gap in seq1
	push(@maplist1,"MISSINGDEN");
    } else {
	$res1one = &one_from_three(substr($reslist1[$j1],7,3));
	if ($res1one ne $seq1[$i]) {
	    print STDERR "Warning: $reslist1[$j1] != $seq1[$i], $j1, $i\n";
	}
	push(@maplist1,$reslist1[$j1++]);
    }

    if ($seq2[$i] eq '-') { # gap in seq2
	push(@maplist2,"MISSINGDEN");
    } else {
	$res2one = &one_from_three(substr($reslist2[$j2],7,3));
	if ($res2one ne $seq2[$i]) {
	    print STDERR "Warning: $reslist2[$j2] != $seq2[$i], $j2, $i\n";
	}
	push(@maplist2,$reslist2[$j2++]);
    }
}


for ($i=0; $i<scalar(@seq1); $i++) {
    print "$maplist1[$i] == $maplist2[$i]\n";
}


exit;
#####








sub one_from_three {
    local $three = shift;
    local $one = ".";
    $one="a" if ($three eq "ALA");
    $one="c" if ($three eq "CYS");
    $one="d" if ($three eq "ASP");
    $one="e" if ($three eq "GLU");
    $one="f" if ($three eq "PHE");
    $one="g" if ($three eq "GLY");
    $one="h" if ($three eq "HIS");
    $one="i" if ($three eq "ILE");
    $one="k" if ($three eq "LYS");
    $one="l" if ($three eq "LEU");
    $one="m" if ($three eq "MET");
    $one="n" if ($three eq "ASN");
    $one="p" if ($three eq "PRO");
    $one="q" if ($three eq "GLN");
    $one="r" if ($three eq "ARG");
    $one="s" if ($three eq "SER");
    $one="t" if ($three eq "THR");
    $one="v" if ($three eq "VAL");
    $one="w" if ($three eq "TRP");
    $one="y" if ($three eq "TYR");
    return uc $one;
}
