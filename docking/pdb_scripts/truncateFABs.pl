#!/usr/bin/perl

# antibody fragments dock at the CDR regions...so the other parts of
# the moledule are not needed.  This script eliminates the unnecessary
# residues: for light chains, residues 112 and up are removed; for
# heavy, 116.

# JJG 1/4/2


# set chain identifiers
$L_ID='L';
$H_ID='H';
if ($ARGV[0] =~ /-(\w)(\w)$/) {
    $chain_arg = shift;
    print STDERR "Chain arg: $chain_arg\n";
    $L_ID=$1;
    $H_ID=$2;
}
elsif ($ARGV[0] =~ /-(\w)$/) {
    $chain_arg = shift;
    print STDERR "Chain arg: $chain_arg  **Dromedary**\n";
    $L_ID='x';
    $H_ID=$1;
}

$L_trunc=112;
$H_trunc=119;

if ($ARGV[0] =~ /-truncL/) {
    shift;
    $L_trunc=shift;
}

if ($ARGV[0] =~ /-truncH/) {
    shift;
    $H_trunc=shift;
}


print STDERR "$0 [-LH|-H] [-truncL x] [-truncH x]\n";
print STDERR "Truncating light chain $L_ID at $L_trunc and heavy chain $H_ID at $H_trunc\n";

print "REMJJG ---Truncated light chain $L_ID at $L_trunc and heavy chain $H_ID at $H_trunc\n";

while (<>) {
    $ch = substr($_,21,1);
    $res_num = substr($_,22,4);
    next if ($ch eq $L_ID && $res_num>=$L_trunc);
    next if ($ch eq $H_ID && $res_num>=$H_trunc); 
    print $_;
}

exit 0;

