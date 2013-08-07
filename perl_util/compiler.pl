#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;

use constant save_fn => '.comp.txt';

my $n_procs = $ENV{N_PROCS};
my $mode    = "release";
my $static  = 0;
my $cxx     = 0;
my $rerun   = 0;
my $app     = "bin pilot_apps_all";

&GetOptions(
	"j=i"    => \$n_procs,
	"m=s"    => \$mode,
	"static" => \$static,
	"c=s"    => \$cxx,
	"a"      => \$app,
	"r"      => \$rerun,
);

my @flags;
if ( $cxx     ) { push @flags, "cxx=$cxx"; }
if ( $n_procs ) { push @flags, "-j$n_procs"; }
if ( $mode    ) { push @flags, "mode=$mode"; }
if ( $static  ) { push @flags, "extras=static"; }
if ( $rerun   ) {
	@flags = @{ read_last_cmd( save_fn ) };
}

my $cmd = join ' ', ('./scons.py',@flags,'bin pilot_apps_all');
print $cmd, "\n";
system( $cmd );

save_last_cmd( save_fn, \@flags );

sub read_last_cmd {
	my $fn = shift;
	open FILE, "<$fn" or die $!;
	my @file = <FILE>;
	close FILE or die $!;

	my @flags = map { chomp $_; $_ } @file;
	return \@flags;
}

sub save_last_cmd {
	my $fn    = shift;
	my $flags = shift;
	open FILE, ">$fn" or die $!;
	foreach my $flag (@flags) {
		print FILE $flag, "\n";
	}
	close FILE or die $!;
}
