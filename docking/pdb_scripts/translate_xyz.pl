#!/usr/bin/perl

# translate_x.pl translates a pdb in the x-direction

# JJG 6/3/2
# modifications made to handle chains of atoms, although can only handle one chain at a time currently.


if ($#ARGV < 4) {
	print STDERR "usage: $0 <deltaX> <deltaY> <deltaZ> <chain> [<pdbfile(s)>]\n";
	print STDERR "move atoms in pdb by a given vector\n";
	exit -1;
}
$dx = shift;
$dy = shift;
$dz = shift;
$chain = shift;

while (<>) {
    if (/^ATOM|^HETATM/) {
	$atom_chain = substr ($_, 21, 1);
	if ($chain eq $atom_chain) {
	    $x = substr ($_, 30, 8);
	    $y = substr ($_, 38, 8);
	    $z = substr ($_, 46, 8);
	    $x += $dx;
	    $y += $dy;
	    $z += $dz;
	    substr($_,30,8) = sprintf ("%8.3f", $x);
	    substr($_,38,8) = sprintf ("%8.3f", $y);
	    substr($_,46,8) = sprintf ("%8.3f", $z);
	}
    }
    print $_;
}

exit 0;
