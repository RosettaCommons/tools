#!/usr/bin/perl

# truncate pdbs at certain residues

# JJG 1/25/2

sub usage() {
    print STDERR "USAGE: $0 chainID gtlt res_num [pdb]\n";
    print STDERR "Removes residues from the specified chain with residues\n";
    print STDERR "greater than (gt) or less than (lt) res_num\n";
    exit;
}

# set chain identifiers
$chainID=shift;
$gtlt=shift;
$trunc=shift || usage();

if ($gtlt eq "gt") {
    $gtltword="above";
} elsif ($gtlt eq "lt") {
    $gtltword="below";
} else {
    print STDERR "gt or lt must be the second argument\n";
    usage();
}
    

print STDERR "Truncating chain $chainID: removing residues $gtltword $trunc\n";

print "REMJJG ---Truncated chain $chainID $gtltword $trunc\n";

while (<>) {
    $card = substr($_,0,4);
    $ch = substr($_,21,1);
    $res_num = substr($_,22,4);
    next if ($card eq "ATOM" && $gtlt eq "gt" && 
	     $ch eq $chainID && $res_num>$trunc);
    next if ($card eq "ATOM" && $gtlt eq "lt" && 
	     $ch eq $chainID && $res_num<$trunc);
    print $_;
}

exit 0;

