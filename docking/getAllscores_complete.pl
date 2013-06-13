#!/usr/bin/perl
# getAllscores_complete.pl: gathers the scorefiles from subdirectories and makes a single scorefile for each pdb ("pdb.sc") and a granddaddy scorefile with all scores in it, labelled with a pdb id ("Allscores")

# adds the pdb in the first column entitled 'id'

if (@ARGV==0) {
# find directories that are four characters and start with a number
    @pdbs = `ls -d ???? | grep "^[1-9|c]"`;
}
else {
    @pdbs = @ARGV;
}

chomp(@pdbs);
warn "Processing: @pdbs\n";
warn "Creating scorefiles/*.sc and scorefiles/Allscores\n";    

# create scorefile directories
`mkdir scorefiles`;
`mkdir backup`;

open(ALLSCORES,">scorefiles/Allscores");

# loop through pdbs
$done=0;
foreach $p (@pdbs) {
    chomp $p; 

# get full atom score files
    @scfiles = `ls $p | grep -E "\.fasc|\.sc" `;

    print "skipping $p\n" if ( @scfiles == 0 );
    next if ( @scfiles == 0 );

# print header once
    if (!$done) { 
	$hdr = `head -1 $p/$scfiles[0]`;
	printf ALLSCORES "id   %s",$hdr;
	$done++;
    }

# open target score file
    open(TARGET,">scorefiles/${p}.sc");
# print target header each time
    $hdr = `head -1 $p/$scfiles[0]`;
    printf TARGET "%s",$hdr;
    @tmp = split(" ",$hdr); $nhdrfields = @tmp;

# loop through all score files in current directory
    foreach $scfi (@scfiles) {
	chomp $scfi;
	printf STDERR "$scfi with $nhdrfields fields"; 

# copy scorefiles for backup
	`cp $p/$scfi backup`;

# filter out non-printable characters
	`cat $p/$scfi | tr -cd "[:print:]|[:space:]" > $p/tmp.sc`;

# open score file for input
	open(II,"${p}/tmp.sc");

	while(<II>) {
	    next if /filename/; # remove header
	    next if /\*\*\*/;   # remove bad entries
	    @tmp = split " "; $nfields = @tmp;
	    next if ($nfields != $nhdrfields); # remove lines with bad columns
	    $frl = $_;
	    printf TARGET "%s",$frl;
## only output the decoys to cout
#	    next unless /output_decoy/;
	    printf ALLSCORES "%4s %s",$p,$frl;
	}
	printf STDERR "\n";
    }
    `rm $p/tmp.sc`;
    close(TARGET);
}

close(ALLSCORES)
#printf STDERR "Remember: cp Allscores scorefiles/\n";
