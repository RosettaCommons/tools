#!/usr/bin/perl -w

use strict;
use File::Find;

# Quick script to add includes to .hh files in the Rosetta source tree. The
# author of this code provides no warranty, no support, and no mercy (implied
# or otherwise).
# James Thompson <tex@u.washington.edu>

my $base_path	 = '/work/tex/src/mini_clean/mini';
#my @source_dirs = qw/ src demo test /;
my @source_dirs = qw/ src /;
#my @source_dirs = qw| src/core |;
my @extensions	= qw/ hh /;

for my $source_dir (@source_dirs) {
	find(
		{ wanted => \&add_include_guards },
		#{ wanted => \&verify_include_guards },
		"$base_path/$source_dir"
	);
}

# makes sure that include guards are generated properly based on the file path
sub verify_include_guards {
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

	my $fullname = $File::Find::dir . "/" . $file;
	# mega-regex hack below!
	$fullname    =~ s/$base_path\/src\///g;
	$fullname    =~ s/fwd\.hh$/FWD\.HH/g;
	$fullname    =~ s/\./\//g;
	$fullname    =~ s/hh$/HH/g;
	my @path     = split /\//, $fullname;
	my $include_guard_symbol = join '_', ('INCLUDED',@path);

	open FILE, "<$file" or die $!;
	my @file = <FILE>;
	close FILE or die $!;

	open FILE, ">$file" or die $!;
	foreach my $line (@file) {
		chomp $line;
		if ( $line =~ /#define (INCLUDED_.*)/ ) {
			#print join ' ', ($fullname,$1,$include_guard_symbol), "\n";
			#if ( $1 ne $include_guard_symbol ) {
			#	print "Error: symbol $1 in $file should be $include_guard_symbol!\n";
			#}
			print FILE "#define $include_guard_symbol\n";
		} else {
			print FILE $line, "\n";
		}
	}
	close FILE or die $!;
}

# Takes a file as an argument and adds include guards to it's #include pragmas.
sub add_include_guards {
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
	open FILE, "<$file" or die $!;
	my @file = <FILE>;
	close FILE or die $!;

	open FILE, ">$file" or die $!;
	my $last_line = '';

	foreach my $line (@file) {
		chomp $line;
		if ( $line =~ /^#include\s+\<([\w\/\.]+)\>/ ) {
			my $include = $1;
			$include =~ s/\./\//g;
			$include =~ s/hh$/HH/g;
			my @path = split /\//, $include;
			my $include_guard_symbol = join '_', ('INCLUDED',@path);

			if ( is_mini_library($path[0]) && $last_line !~ /#ifndef $include_guard_symbol/ ) {
				my $guard_line = "#ifndef $include_guard_symbol"; 
				print FILE $guard_line, "\n"; # include guard
				print FILE $line, "\n";       # original #include pragma
				print FILE "#endif\n";        # end include guard
			} else {
				# already guarded!
				print FILE $line, "\n";
			}
		} else {
			print FILE $line, "\n";
			$last_line = $line;
		}
	}
	close FILE or die $!;
	# end warning.
}

sub is_mini_library {
	my $lib = shift;
	my @valid_libraries = qw/ core protocols numeric utility ObjexxFCL devel /;
	my %is_valid = map { $_ => 1 } @valid_libraries;
	if ( $is_valid{$lib} ) {
		return 1;
	}
	return 0;
}
