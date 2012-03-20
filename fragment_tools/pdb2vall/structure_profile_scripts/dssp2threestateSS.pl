#!/usr/bin/perl -w

use strict;
use File::Basename;
use FindBin qw($Bin);

# Quick program to run DSSP on a PDB file and transform the output
# into three-state secondary structure predictions suitable for
# Rosetta.
# James Thompson <tex@u.washington.edu>

my $dssp_binary = "$Bin/dssp/dssp";

# read config
open(F, "$Bin/../pdb2vall.cfg") or die "ERROR! cannot open config file $Bin/../pdb2vall.cfg: $!\n";
while (<F>) {
  if (/^\s*dssp\s*=\s*(\S+)\s*/) {
    $dssp_binary = $1;
  }
}
close(F);


my $dssp_options  = ' '; #' -na ';

my $usage = <<USAGE
usage: $0 pdbfile
USAGE
;

my $pdb_file = $ARGV[0];
if ( !defined $pdb_file || ! -f $pdb_file ) {
    warn $usage, "\n";
    exit 1;
}

my $dssp_output = `$dssp_binary $dssp_options $pdb_file 2>&1`;
my @output = split /\n/, $dssp_output;

my $outputfile = "$pdb_file.dssp";
open(F, ">$outputfile");
print F $dssp_output;
close(F);

# Correct DSSP output
my @dssp_out;
my $count = 0;

foreach my $line ( @output ) {
    if ($line=~/^\s(....)\s(....)(\s.\s[^!].*)/ and  $count==0) {
        push @dssp_out, " $1 $2$3\n";
    } elsif ($line=~/^\s(....)\s(....)(\s.\s[^!].*)/ and $count!=0) {
        my $formatted = sprintf("%4d",$count);
        push @dssp_out, " $formatted $2$3\n";
        $count++;
    }
    elsif ($line=~/^\s(....)\s(....)(\s.\s!.*)/ and $count==0) {
        print STDERR "CHAIN BREAK $line \n";
        $count=scalar($1);
    } elsif ($line=~/^\s(....)\s(....)\s.\s!.*/ and $count!=0) {
        print STDERR "CHAIN BREAK $line \n";
        #Do nothing
    } elsif ($line!~/^\s(....)\s(....)\s.\s[^!].*/) {
        push @dssp_out, $line;
    } else {
        print $line;
    }
}

# print 3-state secondary structure predictions to STDOUT
print "> DSSP-assigned Secondary Structure for $pdb_file\n";
my $start = 0;
my $i = 0;
foreach my $line (@dssp_out) {
    chomp $line;

    if ($start) {
        my $ss_dssp = substr($line,16,1);
        my $prediction = ss3state($ss_dssp);

        if ( $prediction ne '' ) {
            print $prediction;
            if ( $i % 80 == 79 ) { print "\n" };
            $i++;
        }
    } else {
        if ( $line =~ /\s*#/ ) {
            $start = 1;
        }
    }
}



print "\n";


# ss3state()
# Changes the six-state DSSP alphabet for 2' structure into a
# three-state alphabet. Here's how it works:
# Input         Output
#   H             H
#   I             H
#   G             H
#   E             B
#   B             B
#   L             L
#
# In the shortened alphabet, H denotes a residue that is part of
# an alpha helix, B denotes a residue that is part of a beta sheet,
# and L denotes a residue that is part of loop.

sub ss3state {
    my $ss = shift;

    my $letter;
    if ($ss eq 'H' || $ss eq 'I' || $ss eq 'G') {
        $letter = 'H';
    }
    elsif ($ss eq 'E' || $ss eq 'B') {
        $letter = 'E';
    }
    else {
        $letter = 'L';
    }

    if ( !defined $letter ) {
        warn "Error: can't determine structure for letter ($ss)\n";
    }

    return $letter;
}
