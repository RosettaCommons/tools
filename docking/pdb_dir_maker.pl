#!/usr/bin/perl

# utility to make local directories necessary for rosetta runs.

die("usage: $0 condorscript\n") if ( ! "$ARGV[0]" );

@PDBS = `grep "^pdb=" $ARGV[0]`;

foreach (@PDBS) {
	$pdb_name = substr($_, 4, 4);
	system("mkdir -p $pdb_name/outerr");
}
