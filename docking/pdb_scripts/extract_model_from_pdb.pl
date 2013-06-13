#!/usr/bin/perl

$infile = $ARGV[0];
open(INFILE,"$infile.pdb");
$num=0;
$chain=0;
while ($line=<INFILE>) {
  chomp($line);
  if ($line =~ "MODEL") {
    open(OUTFILE, ">$infile.m.$num.pdb");
    $chain=0;
  } elsif ($line =~ "ENDMDL") {
    print OUTFILE "TER";
    close(OUTFILE);
    $num++;
  } else {
	if (($line =~ " A " && $chain == 2) || ($line =~ " B " && $chain == 1)){
		print OUTFILE "TER\n";
	}	
	if ($line =~ " A "){
		$chain = 1;
	}
	if ($line =~ " B "){
		$chain = 2;
	}
	if ($line !~ "SAM"){
		if ($line !~ "HST"){
			if ($line =~ "ATOM"){
    				print OUTFILE "$line\n";
			}
		}
	}
  }
}
