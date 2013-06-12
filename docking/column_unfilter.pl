#!/usr/bin/perl

# script to check number of columns
die("usage: $0 number_of_columns [files]\n") if (!"$ARGV[0]");
$numCol = shift;

while (<>) {
	chomp $_;
	@w = split /\s+/,$_;
	$numElem = @w;
	if ($numElem != $numCol) {
		print "$_\n";
	}
}
	

