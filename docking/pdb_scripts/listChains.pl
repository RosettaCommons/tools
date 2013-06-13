#!/usr/bin/perl
 
# list chains in the pdb

# JJG 1/3/2


%seen = (); $ch='';

while (<>) {
    if (/^TER/) { push @allchains,"-";}
    next unless (/^ATOM/);
    $ch=substr($_,21,1);
    push @allchains,$ch unless $seen{$ch}++;
}

foreach $ch (@allchains) {
    if ($ch eq "-") { print STDERR "$ch\n"; next; }
    print STDERR "chain $ch $seen{$ch} atoms\n";
}

print @allchains;
