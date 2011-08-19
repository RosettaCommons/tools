#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;


my $verbose = 0;
my $apply_changes = 0;
&GetOptions(
	"apply" => \$apply_changes,
	"verbose" => \$verbose,
);

if ( $apply_changes ) {
	die "You've got to be kidding me.\n";
}

my $fn = 'err';
if (@ARGV) {
	$fn = $ARGV[0];
}

open FILE, "<$fn" or die $!;
my @file = <FILE>;
close FILE or die $!;

my %suggestions_per_file;
foreach my $line (@file) {
	chomp $line;
	if ( $line =~ /(.*)\:(\d+):\s+error:\s+(.*)$/ ) {
		my $file = $1;
		my $line_no = $2;
		my $error_msg = $3;
		if ( $file =~ /owning_ptr/ ) { next; }

		my $suggestion = suggest_solution( $error_msg );
		if ( defined $suggestion ) {
			#print "$file:\n$suggestion\n";
			push @{ $suggestions_per_file{$file} }, $suggestion;
		}
	}
}

foreach my $file (keys %suggestions_per_file) {
	print "$file\n";
	my %printed;
	foreach my $suggestion ( @{ $suggestions_per_file{$file} } ) {
		if ( exists $printed{$suggestion} ) { next; }
		$printed{$suggestion}++;
		print $suggestion, "\n";
	}
	print '-' x 80, "\n";
}

sub suggest_solution {
	my $msg = shift;

	my %func_map = (
		'string_of'     => 'ObjexxFCL/string.functions.hh',
		'I'             => 'ObjexxFCL/fmt/formatted.o.hh',
		'F'             => 'ObjexxFCL/fmt/formatted.o.hh',
		'LJ'            => 'ObjexxFCL/fmt/formatted.o.hh',
		'RJ'            => 'ObjexxFCL/fmt/formatted.o.hh',
		'FArray1D_bool' => 'ObjexxFCL/FArray1D.fwd.hh', 
	);

	my $suggestion;
	if ( $msg =~ /\`(.*)\' was not declared/ ) {
		my $func_name = $1;

		if ( exists $func_map{$func_name} ) {
			my $header = $func_map{$func_name};
			my @tokens = split /\//, $header;
			pop @tokens;
			my $namespace = join '::', @tokens;

			$suggestion .= "#include <$header>\n";
			$suggestion .= "using $namespace" . "::$func_name;\n";
		}
#error: invalid use of undefined type `const struct core::conformation::Residue'
	} elsif ( $msg =~ /invalid use of undefined type \`const struct (.*)'/ ) {
		print "msg = $msg\n";
		my $object_name = $1;
		my $header = $object_name;
		$header =~ s/::/\//g;
		$header .= '.hh';
		$suggestion .= "#include <$header>\n";
	} else {
		if ( $verbose ) {
			print "no suggestion for error: $msg\n";
		}
	}

	return $suggestion;
}
