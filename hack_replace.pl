#!/usr/bin/perl

use strict;
use warnings;

use File::Find;

# hacky script for performing modifications to source files

#usage: edit the script and copy it to rosetta_source/src, then "perl hack_replace.pl".  The mr function lets you specifcy something, but I don't know what.
#to edit the script, first comment on/off the functions you want out of the top make_transformer function.  You'll likely want replace_call and replace_header for XRW-type work.
#next, edit the guts of the functions you'll use to fix the relevant search-and-replace strings.


#In
my $func_ref = make_transformer(
    #\&mr,
				#\&fix_montecarlo_references,
				#\&fix_parsemytag,
				#\&protocol_to_protocols,
				#\&symmetry_fibril_to_fibril,
				#\&prot_symmetry_to_prot_simple_moves
				#\&sfad_to_sd,
				#\&asfd_to_bm,
				#\&basic_moves_to_simple_moves,
				#\&remove_trailing_whitespace,
				\&replace_call,
				\&replace_header,
);

find( { wanted => $func_ref }, '.' );



sub mr {
	my $line = shift;

	$line =~ s#$ARGV[0]#$ARGV[1]#;

	return $line;
}


#this will change header paths, here from old core/io/pdb/file_data.fwd.hh to core/import_pose/file_data.fwd.hh
sub replace_header {
	my $line = shift;

	$line =~ s#<core/io/pdb/file_data.fwd.hh>#<core/import_pose/file_data.fwd.hh>#;

	return $line;
}

#changes function and class names and namespacing
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

	#list the function or class names you want to alter

	my @func_names = qw/
		FileData
	/;
	
	#list all possible namespaces for your function.  Here, core::io::pdb::FileData, so the 4, 3, 2, and 1-unit namespacings need to be considered.
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

sub basic_moves_to_simple_moves {
	my $line = shift;
	$line =~ s/basic_moves/simple_moves/g;
	return $line;
}

sub asfd_to_bm {
  my $line = shift;

  # iguard
  $line =~ s/protocols_moves_asym_fold_and_dock/protocols_basic_moves_asym_fold_and_dock/g;

  #headers
  $line =~ s/protocols\/moves\/asym_fold_and_dock/protocols\/basic_moves\/asym_fold_and_dock/g;

  #namespace
  $line =~ s/protocols::moves::asym_fold_and_dock/protocols::basic_moves::asym_fold_and_dock/g;
  $line =~ s/moves::asym_fold_and_dock/basic_moves::asym_fold_and_dock/g;
  
  return $line;
}

sub sfad_to_sd {
  my $line = shift;

  # iguard
  #$line =~ s/protocols_moves_symmetry_SymFoldandDock/protocols_symmetric_docking_SymFoldAndDock/g;
  #$line =~ s/protocols_moves_symmetry_SetupForSymmetryMover/protocols_symmetric_docking_SetupForSymmetryMover/g;

  #headers
  #$line =~ s/protocols\/moves\/symmetry\/SetupForSymmetryMover/protocols\/symmetric_docking\/SetupForSymmetryMover/g;
  $line =~ s/protocols\/moves\/symmetry\/ExtractAsymmetricUnitMoverCreator/protocols\/symmetric_docking\/ExtractAsymmetricUnitMoverCreator/g;

  #namespace
  #$line =~ s/protocols::moves::symmetry::SymFoldandDock/protocols::simple_moves::symmetry::SymFoldandDock/g;
  #$line =~ s/moves::symmetry::SymFoldandDock/simple_moves::symmetry::SymFoldandDock/g;

  #$line =~ s/protocols::moves::symmetry::SymFoldandDock/protocols::symmetric_docking::SymFoldandDock/g;
  #$line =~ s/moves::symmetry::SymFoldandDock/symmetric_docking::SymFoldandDock/g;

  $line =~ s/protocols::moves::symmetry::ExtractAsymmetricUnitMoverCreator/protocols::symmetric_docking::ExtractAsymmetricUnitMoverCreator/g;
  $line =~ s/moves::symmetry::ExtractAsymmetricUnitMoverCreator/symmetric_docking::ExtractAsymmetricUnitMoverCreator/g;

  return $line;
}

sub prot_symmetry_to_prot_simple_moves {
  my $line = shift;

  # iguard
  $line =~ s/protocols_moves_symmetry_/protocols_simple_moves_symmetry_/g;
  #$line =~ s/protocols_moves_symmetry_SetupForSymmetryMover/protocols_symmetric_docking_SetupForSymmetryMover/g;

  #headers
  #$line =~ s/protocols\/moves\/symmetry\/SetupForSymmetryMover/protocols\/symmetric_docking\/SetupForSymmetryMover/g;
  $line =~ s/protocols\/moves\/symmetry/protocols\/simple_moves\/symmetry/g;

  #namespace
  #$line =~ s/protocols::moves::symmetry::SymFoldandDock/protocols::simple_moves::symmetry::SymFoldandDock/g;
  #$line =~ s/moves::symmetry::SymFoldandDock/simple_moves::symmetry::SymFoldandDock/g;

  #$line =~ s/protocols::moves::symmetry::SymFoldandDock/protocols::symmetric_docking::SymFoldandDock/g;
  #$line =~ s/moves::symmetry::SymFoldandDock/symmetric_docking::SymFoldandDock/g;

  $line =~ s/protocols::moves::symmetry/protocols::simple_moves::symmetry/g;
  $line =~ s/moves::symmetry/simple_moves::symmetry/g;

  return $line;
}

sub symmetry_fibril_to_fibril {
		
  my $line = shift;

  # iguard
  $line =~ s/protocols_moves_symmetry_fibril/protocols_fibril_fibril/g;
  $line =~ s/protocols_moves_symmetry_SetupForFibril/protocols_fibril_SetupForFibril/g;

  #headers
  $line =~ s/protocols\/moves\/symmetry\/fibril/protocols\/fibril\/fibril/g;
  $line =~ s/protocols\/moves\/symmetry\/SetupForFibril/protocols\/fibril\/SetupForFibril/g;

  #namespace
  #$line =~ s/protocols::moves::symmetry::SymFoldandDock/protocols::simple_moves::symmetry::SymFoldandDock/g;
  #$line =~ s/moves::symmetry::SymFoldandDock/simple_moves::symmetry::SymFoldandDock/g;

  $line =~ s/protocols::moves::symmetry::SetupForFibril/protocols::fibril::SetupForFibril/g;
  $line =~ s/protocols::moves::symmetry::fibril/protocols::fibril::fibril/g;

	my @func_names = qw/
		reorient_extended_fibril
		make_symmetric_fibril
		superimpose_pose_on_subset_bb
	/;
	
	foreach my $func_name (@func_names) {
		$line =~ s/protocols::moves::symmetry::($func_name)/protocol::fibril::$1/g;
		$line =~ s/moves::symmetry::($func_name)/protocols::fibril::$1/g;
		$line =~ s/symmetry::($func_name)/protocols::fibril::$1/g;
		$line =~ s/(\s+)($func_name\s*\()/$1protocols::fibril::$2/g;
	}

  return $line;
}

sub protocol_to_protocols {

  my $line = shift;

	my @func_names = qw/
		reorient_extended_fibril
		make_symmetric_fibril
		superimpose_pose_on_subset_bb
	/;

	foreach my $func_name (@func_names) {
		$line =~ s/protocol::fibril::($func_name)/protocols::fibril::$1/g;
	}

  return $line;
}

sub fix_parsemytag {

  my $line = shift;

	$line =~ s/(\s+)DataMap/$1protocols::moves::DataMap/g;
	#$line =~ s/(\s+)Filters_map/$1protocols::moves::Filters_map/g;
	$line =~ s/(\s+)Movers_map/$1protocols::moves::Movers_map/g;
	 
  return $line;
}

sub fix_montecarlo_references{
 my $line = shift;

	#$line =~ s/(\s+)DataMap/$1protocols::moves::DataMap/g;
	#$line =~ s/(\s+)Filters_map/$1protocols::moves::Filters_map/g;
	$line =~ s/(\s+)MonteCarlo/$1protocols::moves::MonteCarlo/g;
	
	$line =~ s/(\s+)RepeatMoverOP/$1protocols::moves::RepeatMoverOP/g;
	
	$line =~ s/^MonteCarlo/protocols::moves::MonteCarlo/g;
	
	$line =~ s/^Mover/protocols::moves::Mover/g;
	 
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
	return ( $_[0] =~ /\.cc$/ || $_[0] =~ /\.hh$/ || $_[0] =~ /IncludeDict/ ||  $_[0] =~ /\.py/ );
}
