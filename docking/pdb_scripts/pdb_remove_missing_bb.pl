#!/usr/bin/perl

# check for the four backbone atoms and remove those with incomplete
# backbones. residue labels must correct, that is, insert letters must
# be used when needed.

# J.Gray 7/02

if ($#ARGV < 0) {
        print STDERR "usage: $0 <pdbfile> \n";
        exit -1;
}

@pdbf=<>;

for $line (@pdbf){

    if ($line =~ /ATOM/) {    
	$reslabel = substr($line,21,6);
	if ($reslabel ne $last_reslabel) { # new residue
	    if ($N && $O && $C && $CA && $CB){
		# complete residue -- add to new pdb
		push(@newpdb,@newatoms);
	    } else { # keep non-ATOM lines always
		for $line2 (@newatoms){
		    push(@newpdb,$line2) if ($line2 !~ /ATOM/);
		}
	    }
	    undef @newatoms;
	    $N=0;
	    $O=0;
	    $C=0;
	    $CA=0;
	    $CB=0;
	    $last_reslabel=$reslabel;
	}
	$atom=substr($line,12,4);
	$three=substr($line,17,3);
	if ($atom eq " N  ") { $N=1 };
	if ($atom eq " O  " ||
	    $atom eq " OT1" ||
	    $atom eq " OT2") { $O=1 };
	if ($atom eq " C  ") { $C=1 };
	if ($atom eq " CA ") { $CA=1 };
	if ($atom eq " CB " ||
	    $three eq "GLY") { $CB=1 };
		    
    }
    push(@newatoms,$line);
}

# last residue
if ($N && $O && $C && $CA && $CB){
    # complete residue -- add to new pdb
    push(@newpdb,@newatoms);
} else {
    for $line2 (@newatoms){
	push(@newpdb,$line2) if ($line2 !~ /ATOM/);
    }
}

print @newpdb;

exit 0;
