#!/usr/bin/perl

# Add insert letters where they are missing in PDBs
# Detect missing insert letters by finding new N atoms
# with the same chain and residue number as the last one
# J.Gray 6/02

if ($#ARGV < 0) {
        print STDERR "usage: $0 <pdbfile> \n";
        exit -1;
}

while (<>){

    $line=$_;
    if (/ATOM/) {
	if (substr($line,12,4) eq ' N  ') { # N = new residue
	    $last_reslabel=$reslabel;
	    $reslabel = substr($line,21,6);
	    if ($reslabel eq $last_reslabel) { # new label required
		$new_reslabel = substr($reslabel,0,5) . $insert_code;
		$insert_code++; # increment the insert letter
		printf STDERR "Relabeling \'$new_reslabel\'\n";
	    }
	}
	if ($reslabel eq $last_reslabel) { # new label required
	    substr($line,21,6) = $new_reslabel;
	} else { # reset the insert letter
	    $insert_code='A';
	}
	    
    }
    print $line;
}
	

exit 0;
