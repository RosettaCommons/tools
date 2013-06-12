#!/usr/bin/perl

#  if (@ARGV < 2) {
#      print STDERR "usage: $0 column cutoff datafile\n";
#      exit;
#  }

$column=shift;
$cutoff=shift || die("usage: $0 column cutoff datafile\n");

while (<>) {

    @fields = split " ";

    if (/filename/) {
	# header line
	printf STDERR "filtering $datafile on $fields[$column] < $cutoff\n";
	printf $_;
	next;
    }

#    print $fields[0]," ",$fields[5],"\n";

    next if ($fields[$column] > $cutoff);

    printf $_;

}

#print "\n";
