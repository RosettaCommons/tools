#!/usr/bin/perl

if ($#ARGV < 0) {
        print STDERR "usage: $0 <pdbfile> \n";
        exit -1;
}

$pdbfile=$ARGV[0];

@residues=`grep ^ATOM $pdbfile | cut -c 18-27 |uniq`;
@atomcounts=`grep ^ATOM $pdbfile | cut -c 18-27 |uniq -c |cut -f1`;
chomp(@residues);

$natoms=0;
$nres=0;
$natoms_total=0;
$nres_total=0;
$i=0;
for $r (@residues) {
    $ch = substr($r,4,1);
    $n = substr($r,5,4);
#    print "ch $ch n $n --$r";

    if ($ch ne $last_ch || $n != $last_n + 1) {
	print "$first_r -- $last_r: $nres residues; $natoms atoms\n" if ($first_r);
	$first_r = $r;
	chomp($first_r);
	$natoms=0;
	$nres=0;
    }

    $last_ch = $ch;
    $last_n  = $n;
    $last_r  = $r;
    $nres++;
    $nres_total++;
    $natoms += $atomcounts[$i];
    $natoms_total += $atomcounts[$i++];
}

print "$first_r -- $last_r: $nres residues; $natoms atoms\n";
print "                  TOTAL : $nres_total residues; $natoms_total atoms\n";
