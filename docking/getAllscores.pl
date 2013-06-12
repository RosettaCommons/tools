#!/usr/bin/perl
# this is to go through and get all scores to do plot scores
# also does reweight according to logistic regression of 12.04 runs
# needs pdb ss type from runlist

# adds the pdb in the first column entitled 'id'

if (@ARGV==0) {
# find directories that are four characters and start with a number
    @pdbs = `ls -d ???? ????? | grep "^[1-9|c]"`;
}
else {
    @pdbs = @ARGV;
}

warn "Processing: @pdbs\n";    

# create scorefile directories
`mkdir backup`;

# loop through pdbs
$done=0;
foreach $p (@pdbs) {
    chomp $p; 

# get full atom score files
    @scfiles = `ls $p | grep "\.[fa|]sc" `;

# print header once
    if (!$done) { 
	$hdr = `head -1 $p/$scfiles[0]`;
	printf "id   %s",$hdr;
	$done++;
    }

# loop through all score files in current directory
    foreach $scfi (@scfiles) {
	chomp $scfi;
	printf STDERR "$scfi";

# copy scorefiles for backup
	`cp $p/$scfi backup`;

# open score file for input
	open(II,"${p}/$scfi");

	while(<II>) {
	    next if /\*\*\*/;  # remove bad entries
	    next unless /output_decoy/; # only output the decoys to cout
	    printf "%4s %s",$p,$_;
	}
	printf STDERR "\n";
    }
}
