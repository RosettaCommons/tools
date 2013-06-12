#!/usr/bin/perl

sub usage() {
    print "USAGE: $0 column_number ge/le limit [file]\n";
    exit;
}

$column = shift;
$gele = shift;
$limit = shift || usage();

if ($gele eq "ge") {
    $geleword="equal or more than";
} elsif ($gele eq "le") {
    $geleword="equal or less than";
} else {
    print STDERR "ge or le must be the second argument\n";
    usage();
}
    
print STDERR "Finding the rank of the first decoy with column $column $geleword $limit\n";

$header=<>;
$i=0;
while (<>) {
    @fields=split(/\s+/,$_); # split at whitespace
    next if (@fields < $column);
    $i++;
    $x = $fields[$column]; #print STDERR "$x $limit\n";
    if (($gele eq "ge" && $x >= $limit) ||
	($gele eq "le" && $x <= $limit)) {
	print "\t$i";
	print STDERR "decoy $i has $x $geleword $limit\n";
	exit;
    }
}

print "\t>$i";
exit;


###
#findRank.pl 17 ge 0.25 clustersummary  ## examines contact ratio
#findRank.pl 5 le 5 clustersummary      ## examines rms
