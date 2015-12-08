#!/usr/bin/perl

# find the rmsd between all pairs of decoys in the current directory
# usage: rms2.pl < list_of_pdbs > results_tabular
# jjg 9/5/1

# rev. 9/14: now tolerant of missing/broken pdb files/links


# get list of pdbs
#@pdbs = `find ./ -name "*pdb"`;
#chomp(@pdbs);
#print STDOUT "argv[1]: $ARGV[0]\n";

if (@ARGV < 1) {
    print STDERR "Usage: $0   list_of_pdbs > results_tabular  OR\n";
    print STDERR "Usage: $0 < list_of_pdbs > results_tabular\n";
    @pdbs = <STDIN>;
}
else {
    open (IN,$ARGV[0]);
    @pdbs = <IN>;
#	print "pdbs: @pdbs\n";
};

chomp(@pdbs);

# store CA locations
for ($i=0; $i < @pdbs; $i++) {
    ($ca[$i]) = [&getcoords($pdbs[$i])];
	$coors = &getcoords($pdbs[$i]);
    if (!$ca[$i][0]) {
	# remove pdbs that do not have coordinates
#	print STDERR "CAUGHT!\n";
	splice(@pdbs,$i,1);
	$i--;
	next;
    }
    $rms[$i][$i]=0;
}


# print header
print "pdbs ";
foreach $pdb (@pdbs) { print $pdb," "; }
print "\n"; #,scalar(@pdbs),"\n";


# do rmsds
for ($i=0; $i < @pdbs-1; $i++) {
    for ($j=$i+1; $j < @pdbs; $j++) {
	printf STDERR "%s(%i)-%s(%i) ",$pdbs[$i],$i,$pdbs[$j],$j;
	$rms[$i][$j] = &rms($ca[$i],$ca[$j]);
	$rms[$j][$i] = $rms[$i][$j];
	printf STDERR "rmsd: $rms[$i][$j]\n";
    }
}

# print matrix
for ($i=0; $i < @pdbs; $i++) {
    print $pdbs[$i]," ";
    for ($j=0; $j < @pdbs; $j++) {
	printf "%10.4f ",$rms[$i][$j];
    }
    print "\n";
}

# done!


#print @pdbs;

#($ca[1]) = [&getcoords("oDcpr2s.ppk_0015.pdb")];
#($ca[2]) = [&getcoords("obcpr2s.ppk_0010.pdb")];

#print $ca[1];

#print rms($ca[1],$ca[0]),"\n";

#print $ca[1][4],"\n";



#----------------------------------------------------
# get coordinates from the given pdb file 
# (CA only, second docking partner)
sub getcoords {
    my $pdbfile = @_[0];
    my (@coords,$part2,@fields);
    if (!open(PDB,$pdbfile)){
	print STDERR "file $pdbfile not found\n";
	@retval = 0;
	return;
    }
    
#    print STDERR "$pdbfile\n";
    $part2=0;
    while(<PDB>) {
		#next unless (/TER/ || $part2);
		#$part2 = 1; # found part2
	next unless /CA/;
	@fields = split ' ';
#	push(@coords,$fields[6],$fields[7],$fields[8]);
	push(@coords,substr($_,30,8),substr($_,38,8),substr($_,46,8));
#	print "fields678: ",$fields[6],$fields[7],$fields[8],"\n";
#	print "substr   : ",substr($_,30,8),substr($_,38,8),substr($_,46,8),"\n";
    }
    
    $clen = @coords;
#    print "coords: ",$coords[4]," ",$clen,"\n";
    @retval = @coords;
}
    


#----------------------------------------------------
# rms between two sets of inputted coordinates
sub rms {
    my ($coords1,$coords2) = @_;
    my ($i,$sum,$rmsd,$retval);
    
    $sum = 0;
    for ($i=0; $i < @$coords1; $i++) {
	$sum = $sum + ($$coords1[$i] - $$coords2[$i])**2;
    }

    $rmsd = sqrt($sum / (@$coords1/3) );
#    print STDERR "rmsd: ",$rmsd,"\n";
    $retval = $rmsd;
}
