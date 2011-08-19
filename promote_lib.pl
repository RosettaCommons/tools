#!/usr/bin/perl

use strict;
use warnings;

use Cwd qw/ abs_path /;
use File::Find;
use File::Path qw/ mkpath /;
use Getopt::Long;

# configuration
my %options;
$options{src_path}   = "/work/tex/src/mini/src";
$options{orig_lib}   = "devel";
$options{dest_lib}   = "protocols";
$options{lib_name}   = "NoesyAssign";
$options{extensions} = [ 'hh', 'cc' ];

&GetOptions(
	\%options,
	"src_path=s",
	"orig_lib=s",
	"dest_lib=s",
	"lib_name=s",
	"extensions=s",
);

my $orig_dir = assemble_path( $options{src_path}, $options{orig_lib}, $options{lib_name} );
my $dest_dir = assemble_path( $options{src_path}, $options{dest_lib}, $options{lib_name} );

if ( ! -d $orig_dir ) {
	die "Error: directory $orig_dir doesn't exist!\n";
}

# maybe this isn't necessary ...
if ( -d $dest_dir ) {
	die "Error: directory $dest_dir already exists!\n";
}
mkpath $dest_dir;

my $func_ref = make_renamer( \%options );
find( $func_ref, $options{src_path} );

sub make_renamer {
	my $options = shift;

	my $src_path   = $options->{src_path};
	my $orig_lib   = $options->{orig_lib};
	my $dest_lib   = $options->{dest_lib};
	my $lib_name   = $options->{lib_name};
	my @extensions = @{ $options->{extensions} };
	
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

		#$file = $File::Find::dir . '/' . $file;
		my @new_lines;
		my $n_changes = 0;
		LINE: foreach my $line (@file) {
			chomp $line;
			if ( !defined $line ) {
				die "Error: line $line not defined in $file!\n";
			}

			if ( $line =~ /^\s*$/ ) {
				push @new_lines, $line;
			} elsif ( $line =~ /^\/\// ) {
				# skip comments
				push @new_lines, $line;
			} elsif ( $line =~ /namespace\s*(\w+)/ && $1 eq $orig_lib && $File::Find::dir =~ /$lib_name/ ) {
				# rename namespace
				my $copy = $line;
				$copy =~ s/$orig_lib/$dest_lib/;
				push @new_lines, $copy;
				$n_changes++;
			} elsif ( $line =~ /\#include/ && $line =~ /$orig_lib\/$lib_name/ ) {
				# rename #include directives
				my $copy = $line;
				$copy =~ s/$orig_lib/$dest_lib/;
				push @new_lines, $copy;
				$n_changes++;
			} elsif ( $line =~ /^(.*)\s+(INCLUDED_.*$orig_lib\_$lib_name\_.*)/ ) {
				# rename #include guards
				my $directive = $1;
				my $tag = $2;
				my $new_tag = $tag;
				$new_tag =~ s/$orig_lib/$dest_lib/g;
				#my $new_line = join ' ', ( $directive, $new_tag );
				push @new_lines, join ' ' , ( $directive, $new_tag );	
				$n_changes++;
			} elsif ( $line =~ /$orig_lib\:\:$lib_name/ ) {
				my $new_line = $line;
				# replace function calls
				$new_line =~ s/$orig_lib\:\:$lib_name/$dest_lib\:\:$lib_name/g;
				push @new_lines, $new_line;
				$n_changes++;
			} else {
				push @new_lines, $line;
			}
		}

		my $dest_fn = assemble_path($File::Find::dir,$file);
		if ( abs_path($File::Find::dir) eq abs_path(assemble_path($src_path,$orig_lib,$lib_name)) ) {
			$dest_fn = assemble_path( $src_path, $dest_lib, $lib_name, $file );
		}

		#if ( -f $dest_fn ) {
		#	die "Error: not overwriting $dest_fn as it already exists!\n";
		#}
		if ( $n_changes > 0 ) {
			print "writing $dest_fn ... ";
			open FILE, ">$dest_fn" or die $!;
			foreach my $line (@new_lines) {
				print FILE $line, "\n";
			}
			close FILE or die $!;
			print "done.\n";
		}
	};
}

sub assemble_path {
	return join '/', @_;
}
