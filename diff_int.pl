#!/usr/bin/perl

use strict;
use warnings;

my $output = `./integration.py --compareonly`;

my @lines = split /\n/, $output;

foreach my $line (@lines) {
	chomp $line;
	if ( $line =~ /^\s*Files (.*) and (.*) differ$/ ) {
		#print $line, "\n";
		my $ref_fn = $1;
		my $new_fn = $2;
		my $output = `diff $1 $2`;
		print $output, "\n";
	}
}
