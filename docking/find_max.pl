#!/usr/bin/perl

$datafile=$ARGV[0];
$column=$ARGV[1];
$cutoff=$ARGV[2];

if (@ARGV < 2) {
    print STDERR "usage: $0 datafile column\n";
    exit;
}

open(SCFILE,$datafile);

$maximum = -99999.9;
while (<SCFILE>) {

    @fields = split " ";

    if (/score/) {
	# header line
	printf STDERR "searching $datafile for maximum of $fields[$column]...";
	next;
    }

#    print $fields[0]," ",$fields[5],"\n";

    next if ($fields[$column] < $maximum);
    next if ($fields[$column] eq "");

    $maximum = $fields[$column];

}

printf STDERR "$maximum\n";
printf "$maximum\n";

