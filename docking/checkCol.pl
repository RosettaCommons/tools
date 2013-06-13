#!/usr/bin/perl

# script to check number of columns
$file = $ARGV[0];
$numCol = $ARGV[1];
open(INFILE,"$file") || die " couldn't find $file\n";
while( $line = <INFILE>) {
	chomp $line;
	@w = split /\s+/,$line;
#	unless ($w[0]) {$numCol++;} ## split will make $w[0] empty is leading whitespace on line
	$numElem = @w;
	if ($numElem == $numCol) {
		print "$line\n";
	}
}
	

