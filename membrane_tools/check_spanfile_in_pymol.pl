#!/usr/bin/env perl

#(c) Copyright Rosetta Commons Member Institutions.
#(c) This file is part of the Rosetta software suite and is made available under license.
#(c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
#(c) For more information, see http://www.rosettacommons.org. Questions about this can be
#(c) addressed to University of Washington CoMotion, email: license@uw.edu.

# author: JKLeman (julia.koehler1982@gmail.com)

$spanfile = $ARGV[0]; 	#input
$pdbfile = $ARGV[1]; 	#input
$out = $spanfile.".pml"; #output
$pdb = substr($pdbfile, 0, 4);

$nres = `grep CA $pdbfile | wc -l`;
$color_sol = "blue";
$color_tm = "orange";

if ($#ARGV < 1){
	print "usage: check_spanfile_in_pymol.pl <spanfile> <pdbfile>\n";
}

#print header for pymol
open ( IN, $spanfile ) || die "cannot open ".$spanfile."!\n";
open ( OUT, ">".$out ) || die "cannot open output file".$pdbfile."!\n";
print OUT "delete everything\n";
print OUT "bg white\n";
print OUT "set ray_shadows=0\n";
print OUT "load ".$pdbfile.", ".$pdb."\n";
print OUT "hide everything\n";
print OUT "color ".$color_sol.", ".$pdb."\n";
print OUT "show cartoon, ".$pdb."\n";

#put spans from spanfile in arrays
$lc = 0;
while ( $in = <IN> ){
	chomp $in;
	$lc++;
	
	# get rid of multiple whitespaces
	$in =~ s/\s+//;

	# split at whitespace
	@in = split /\s+/, $in;

	# original format
	if ( $lc > 4 && $in !~ /^[A-Za-z_]/ ) {

		print "original format\n";
		
		# get start / end residues
		push @start, $in[0];
		push @end, $in[1];
	} 
	# format with chain and insertion code
	elsif ( $lc > 4 && $in =~ /^[A-Za-z_]/ ) {

		print "format with chain and insertion code\n";

		# get chain: first character in string
		push @chain, substr( $in[0], 0, 1 );
		
		# get rid of insertion code
		$icode1 = chop( $in[0] );
		$icode2 = chop( $in[1] );
		
		# get start / end residues
		push @start, substr( $in[0], 1, 10 );
		push @end, substr( $in[1], 1, 10 );
	}
}
close IN;

#go through residues, then check with arrays
$i = 0;
foreach $r ( 1..$nres ){
#	print "$r\t$i\t$start[$i]\t$end[$i]\n";

	if ( $r >= $start[$i] and $r <= $end[$i] ){
		
		# original format
		if ( scalar(@chain) == 0 ) {
			print OUT "color ".$color_tm.", resi $r\n";
		}
		# format with chain
		if ( scalar(@chain) > 0 ) {
			print OUT "color ".$color_tm.", resi $r and chain $chain[$i]\n";
		}

		# increment spans
		if ( $r == $end[$i] ){
			$i++;
		}
	}
}
close OUT;
print "pymol $out &\n";
system ("pymol $out &");
