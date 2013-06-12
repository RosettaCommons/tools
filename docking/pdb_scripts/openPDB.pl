#!/usr/bin/perl

# openPDB: separate the two docking components

# JJG 10/19/01



if ($#ARGV < 1) {
	print STDERR "usage: $0 <pdbfile> <distance> \n";
	print STDERR "move docking partner 2 by a given distance\n";
	print STDERR "along the line of centers to open up the complex\n";
	exit -1;
}
$pdbfile = shift @ARGV;
$distance = shift @ARGV;

open (PDB, $pdbfile);

# first loop: find centroid of partner 1
while (<PDB>) {
    if (/^TER/) {
	push (@buf, $_);
	last;
    }
    if (! /^ATOM/) {
	push (@buf, $_);
	next;
    }
    $atom = substr ($_, 13, 3);
    if ($atom !~ /^H/) {
	$x = substr ($_, 30, 8);
	$y = substr ($_, 38, 8);
	$z = substr ($_, 46, 8);
	$x_sum += $x;
	$y_sum += $y;
	$z_sum += $z;
	++$atoms;
    }  
    push (@buf, $_);
}

$x_av1 = $x_sum / $atoms;
$y_av1 = $y_sum / $atoms;
$z_av1 = $z_sum / $atoms;

$partner2_start = @buf; #print STDERR "2start: $partner2_start\n";
$atoms=0;
$x_sum=0.0;
$y_sum=0.0;
$z_sum=0.0;

# second loop: find centroid of partner 2
while (<PDB>) {
    if (/^TER/) {
	push (@buf, $_);
	last;
    }
    if (! /^ATOM/) {
	push (@buf, $_);
	next;
    }
    $atom = substr ($_, 13, 3);
    if ($atom !~ /^H/) {
	$x = substr ($_, 30, 8);
	$y = substr ($_, 38, 8);
	$z = substr ($_, 46, 8);
	$x_sum += $x;
	$y_sum += $y;
	$z_sum += $z;
	++$atoms;
    }  
    push (@buf, $_);
}
close (PDB);

$x_av2 = $x_sum / $atoms;
$y_av2 = $y_sum / $atoms;
$z_av2 = $z_sum / $atoms;

# calculate displacement vector
$x_loc = $x_av2-$x_av1; # line of centers
$y_loc = $y_av2-$y_av1;
$z_loc = $z_av2-$z_av1;

$d_loc = sqrt($x_loc**2+$y_loc**2+$z_loc**2);

$x_delta = $x_loc * $distance / $d_loc;
$y_delta = $y_loc * $distance / $d_loc;
$z_delta = $z_loc * $distance / $d_loc;

for ($i=$partner2_start; $i <= $#buf; ++$i) {
	next if ($buf[$i] !~ /^ATOM/);

	$x = substr ($buf[$i], 30, 8);
	$y = substr ($buf[$i], 38, 8);
	$z = substr ($buf[$i], 46, 8);

	$x = $x + $x_delta;
	$y = $y + $y_delta;
	$z = $z + $z_delta;  
	
	substr ($buf[$i], 30, 8) = sprintf ("%8.3f", $x);
	substr ($buf[$i], 38, 8) = sprintf ("%8.3f", $y);
	substr ($buf[$i], 46, 8) = sprintf ("%8.3f", $z);
}


print @buf;
exit 0;

#######
# end #
#######
