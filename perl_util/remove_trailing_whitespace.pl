#!/usr/bin/perl -w

use strict;
use File::Find;
use Getopt::Long;

my $quiet = 0;
&GetOptions(
	"quiet" => \$quiet
);

# Quick script to remove trailing whitespace from Rosetta source tree. The
# author of this code provides no warranty, no support, and no mercy (implied
# or otherwise).
# James Thompson <tex@u.washington.edu>

my $base_path;
if ( @ARGV ) {
	$base_path = $ARGV[0];
} else {
	#backwards compatable for hopefully more people than james
	$base_path = $ENV{HOME} . '/src/mini';
	if( not -e $base_path ) {
		$base_path = '.';
	}
}

my @source_dirs = qw/ src demo test tools /;
my @extensions	= qw/ cc hh py /;

my @problem_regexes = (
	'std::cout',
	'std::cerr'
);

my @problem_lines;
for my $source_dir (@source_dirs) {
	find(
		{ wanted => \&remove_trailing_whitespace },
		"$base_path/$source_dir"
	);
}

print "problems with source code that someone should fix:\n";
foreach my $line (@problem_lines) {
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

	# warning: changing the following code will void your warranty!
	if ( !$quiet ) {
		print "removing whitespace from file $file ... ";
	}
	open FILE, "<$file" or die $!;
	my @file = <FILE>;
	close FILE or die $!;

	open FILE, ">$file" or die $!;
	my $count = 1;
	foreach my $line (@file) {
		chomp $line;
		$line =~ s/\s+$//;
		print FILE $line, "\n";

		foreach my $regex (@problem_regexes) {
			if ( $line =~ $regex ) {
				# check for comments
				if ( $line !~ /\s*\/\//) {
					push @problem_lines, "$file:$count:$line";
				}
			}
		}
		$count++;
	}
	close FILE or die $!;
	if ( !$quiet ) {
		print "done.\n";
	}
	# end warning.
}
