#!/usr/bin/perl -w

use strict;
use File::Find;
use Getopt::Long;
use File::Basename;

$0 = basename $0;

my $help = 0;
my $analyze = 0;
my $base_path = $ENV{PWD};
&GetOptions(
	'path=s'  => \$base_path,
	'analyze' => \$analyze,
	'help'    => \$help,
);

my $usage = <<USAGE;
usage:  $0 [options]
  options are:
	--path <path_to_source>
	--analyze
	--help
USAGE

if ( ! -d $base_path ) {
	print STDERR "Error: $base_path not defined!\n";
	print STDERR $usage, "\n";
	exit 1;
}

if ( $help ) {
	print STDERR $usage, "\n";
	exit 0;
}

my %include_map;
my $func_ref = make_func_ref( \%include_map );
my @extensions	= qw/ cc hh /;
find(
	{ wanted => $func_ref },
	$base_path
);

if ( $analyze ) {
	my %ratio;
	print join ' ', qw/ file includes fwd_includes ratio /, "\n";
	foreach my $file (keys %include_map) {
		if ( $file !~ /\.hh$/ ) { next; }
		my $n_includes = scalar @{ $include_map{$file} };
		my $n_fwd_includes = scalar
			grep { $_ =~ /\.fwd\.hh/ }
			@{ $include_map{$file} };

		my $ratio = $n_fwd_includes / $n_includes;
		print join ' ', (
			$file,
			$n_includes,
			$n_fwd_includes,
			$ratio,
		);
		print "\n";
	}
} else {
	foreach my $file (keys %include_map) {
		foreach my $included_file ( @{ $include_map{$file} } ) {
			print join ': ', ($file, $included_file);
			print "\n";
		}
	}
}

sub make_func_ref {
	my $map = shift;

	my $func_ref = sub {
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
		my $skip_regex = qr{ #ifdef WIN32 };
		my $skip_count = 0;
		foreach my $line (@file) {
			if ( $line =~ /\/\// ) { next; } # skip comments
			if ( $line =~ $skip_regex ) {
				$skip_count = 1;
			}
			if ( $line =~ /#endif/ && $skip_count >= 1 ) {
				$skip_count--;
			}
			if ( $line =~ /#include/ && $skip_count == 0 ) {
				if ( $line =~ /^\s*#include\s+<([\d\w\/\.]+)>/ ) {
					#print "$file: $line";
					push @{ $map->{$file} }, $1;
				} elsif ( $line =~ /#include\s+"([\d\w\/\.]+)"/ ) {
					push @{ $map->{$file} }, $1;
				} else {
					warn "Don't recognize include on $file:$line!\n";
					push @{ $map->{$file} }, $line;
				}
			}
		}
	};

	return $func_ref;
}
