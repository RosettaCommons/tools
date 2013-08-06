#!/usr/bin/perl -w

use strict;
use File::Find;
use Getopt::Long;

&GetOptions();

my $base_path	 = '/home/andrew/rosetta/SVN/mini';
my @source_dirs = qw/ src demo test tools /;
my @extensions	= qw/ cc hh py /;

my @strings = (
	'SingleResidueFragData'
);
my @regexes = map { qr{$_} } @strings;

#my $str = "#include <core/scoring/ScoreFunction.cc>";
#foreach my $r (@regexes) {
#	if ( $str =~ $r ) {
#		print "$str matches $r!\n";
#	} else {
#		print "$str doesn't match $r!\n";
#	}
#}

my @found_lines;
for my $source_dir (@source_dirs) {
	find(
		{ wanted => \&remove_trailing_whitespace },
		"$base_path/$source_dir"
	);
}

foreach my $line (@found_lines) {
	print $line, "\n";
}


# Takes a file as an argument and removes all trailing whitespace from its
# lines.
sub remove_trailing_whitespace {
	my $file = $_;

	if ( $file =~ /\.svn/ ) {
		#print "ignoring $file\n";
		return 0;
	}

	if ( ! -f $file ) {
		#print "$file isn't a file!\n";
		return 0;
	}

	my $success = 0;
	for my $extension ( @extensions ) {
		if ( $file =~ /\.$extension$/ ) {
			$success = 1;
			last;
		}
	}
	if ( ! $success ) {
		return 0;
	}

	open FILE, "<$file" or die $!;
	my @file = <FILE>;
	close FILE or die $!;

	$file = $File::Find::dir . '/' . $file;
	my $count = 1;
	foreach my $line (@file) {
		chomp $line;

		foreach my $regex (@regexes) {
			if ( $line =~ $regex ) {
				# check for comments
				if ( $line !~ /\s*\/\//) {
					push @found_lines, "$file:$count:$line";
				}
			}
		}
		$count++;
	}
	# end warning.
}
