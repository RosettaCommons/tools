#!/usr/bin/perl

use strict;
use warnings;

use File::Find;

# hacky script for performing modifications to source files

my $func_ref = make_transformer(
	\&remove_trailing_whitespace,
	 \&replace_call,
	#\&replace_header,
);

find( { wanted => $func_ref }, '.' );


sub replace_header {
	my $line = shift;

	$line =~ s#<core/io/pdb/file_data.fwd.hh>#<core/import_pose/file_data.fwd.hh>#;

	return $line;
}

sub replace_call {
	my $line = shift;
	#my @func_names = qw/
		#centroid_pose_from_pdb
		#poseOPs_from_pdbs
		#pose_from_pdb
		#pose_from_pdbstring
		#pose_from_pose
		#poses_from_pdbs
		#read_additional_pdb_data
		#write_additional_pdb_data
	#/;
	my @func_names = qw/
		FileData
	/;
	
	foreach my $func_name (@func_names) {
		$line =~ s/core::io::pdb::($func_name)/core::import_pose::$1/g;
		$line =~ s/io::pdb::($func_name)/core::import_pose::$1/g;
		$line =~ s/pdb::($func_name)/core::import_pose::$1/g;
		$line =~ s/(\s+)($func_name\s*\()/$1core::import_pose::$2/g;
	}

	return $line;
}

sub remove_trailing_whitespace {
	my $line = shift;
	$line =~ s/\s+$//g;
	return $line;
}

sub basic_moves_to_basic_moves {
	my $line = shift;
	$line =~ s/basic_moves/basic_moves/g;
	return $line;
}


sub make_transformer {
	my @funcs = @_;

	return sub {
		my $fn = $_;
		#print "reading $fn\n";

		if ( is_source_file($fn) ) {
			# read file
			open FILE, "<$fn" or die $!;
			my @lines = <FILE>;
			close FILE or die $!;

			# write file
			my $backup = 0;
			open FILE, ">$fn" or die $!;
			foreach my $line (@lines) {
				my $copy = $line;
				chomp $copy;
				foreach my $func (@funcs) {
					$copy = &$func($copy);
				}

				if ( $copy ne $line ) {
					$backup = 1;
				}
				print FILE $copy, "\n";;
			}
			close FILE or die $!;
			#die "Finished with $fn\n";

			# backup if necessary
			#if ( $backup ) {
			#	open FILE, ">$fn.bak" or die $!;
			#	foreach my $line (@lines) {
			#		print FILE $line;
			#	}
			#	close FILE or die $!;
			#}
		} # is_source_file
	}
}

sub is_source_file {
	return ( $_[0] =~ /\.cc$/ || $_[0] =~ /\.hh$/ );
}
