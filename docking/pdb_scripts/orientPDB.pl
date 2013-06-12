#!/usr/bin/perl

# orientPDB: given a residue, orient the PDB so that the vector from
# the protein center to that residue is oriented parallel to the
# positive x-axis

# built on dylan's rotatePDB script

# JJG 7/27



if ($#ARGV < 2) {
	print STDERR "usage: $0 <pdbfile> <res> <chain> <[center] [flip]>\n";
	print STDERR "rotate the given residue to point along the x-axis\n";
	print STDERR "center will move the centroid to the x-axis\n";
	print STDERR "flip will orient the residue the other way\n";
	exit -1;
}
$pdbfile = shift @ARGV;
$res = shift @ARGV;
$chain = shift @ARGV;
$center = shift @ARGV;
$flip = shift @ARGV;

$degtorad = 3.14159265358 / 180.0;
$radtodeg = 180.0 / 3.14159265358;

open (PDB, $pdbfile);
while (<PDB>) {
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

	if (substr ($_, 22, 4) == $res && 
	    $atom =~ /CA/ && 
	    substr($_, 21, 1) eq $chain ) {
	    print STDERR $_;
	    $found_res = 1;
	    $x0 = $x;
	    $y0 = $y;
	    $z0 = $z;
	}
	    
	push (@buf, $_);
}
close (PDB);

if (!$found_res) {
    print STDERR "Residue ",$res," not found\n";
    exit;
}


$x_av = $x_sum / $atoms;
$y_av = $y_sum / $atoms;
$z_av = $z_sum / $atoms;

$x0 = $x0 - $x_av;
$y0 = $y0 - $y_av;
$z0 = $z0 - $z_av;

if (!$flip) {
    print STDERR "flipping\n";
    $x0 = - $x0;
    $y0 = - $y0;
    $z0 = - $z0;
}

$zrot = - atan2 ($y0,$x0);
$yrot =   atan2 ($z0,sqrt($x0*$x0+$y0*$y0));
		
print STDERR $x0," ",$y0," ",$z0,"\n";
print STDERR "rotations-- z:",$zrot*$radtodeg," y:",$yrot*$radtodeg,"\n";

$zrotmat = [[  cos $zrot, -sin $zrot,  0 ],
	    [  sin $zrot,  cos $zrot,  0 ],
	    [         0,         0,  1 ]];

$yrotmat = [[  cos $yrot,  0,  sin $yrot ],
	    [         0,  1,         0 ],
	    [ -sin $yrot,  0,  cos $yrot ]]; 


for ($i=0; $i <= $#buf; ++$i) {
	next if ($buf[$i] !~ /^ATOM/);

	$x = substr ($buf[$i], 30, 8);
	$y = substr ($buf[$i], 38, 8);
	$z = substr ($buf[$i], 46, 8);

	$x_rel = $x - $x_av;
	$y_rel = $y - $y_av;
	$z_rel = $z - $z_av;

	($x_rel, $y_rel, $z_rel) = &rotate ($zrotmat, $x_rel, $y_rel, $z_rel);
	($x_rel, $y_rel, $z_rel) = &rotate ($yrotmat, $x_rel, $y_rel, $z_rel);

	if ($center) {
	    $x = $x_rel + $x_av;
	    $y = $y_rel;
	    $z = $z_rel;
	}
	else
	{
	    $x = $x_rel + $x_av;
	    $y = $y_rel + $y_av;
	    $z = $z_rel + $z_av;
	}
	    

	substr ($buf[$i], 30, 8) = sprintf ("%8.3f", $x);
	substr ($buf[$i], 38, 8) = sprintf ("%8.3f", $y);
	substr ($buf[$i], 46, 8) = sprintf ("%8.3f", $z);
}

# translate the axes (testing...)

if (!$center) {
    for ($i=0; $i <= $#buf; ++$i) {
	next if ($buf[$i] !~ /^HETATM/);
	
	$x = substr ($buf[$i], 30, 8);
	$y = substr ($buf[$i], 38, 8);
	$z = substr ($buf[$i], 46, 8);
	
	$y = $y + $y_av;
	$z = $z + $z_av;
	
       	substr ($buf[$i], 30, 8) = sprintf ("%8.3f", $x);
	substr ($buf[$i], 38, 8) = sprintf ("%8.3f", $y);
	substr ($buf[$i], 46, 8) = sprintf ("%8.3f", $z);
    }
}


print @buf;
exit 0;

########
# subs #
########

sub rotate {
	local ($rotmat, @cart) = @_;
	local @out = (0, 0, 0);
	my ($i, $j);
	for ($i=0; $i < 3; ++$i) {
		for ($j=0; $j < 3; ++$j) {
			$out[$i] += $rotmat->[$i]->[$j] * $cart[$j];
		}
	}
	return @out;
}

#######
# end #
#######
