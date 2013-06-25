#!/usr/bin/perl
use FindBin qw( $Bin );
$| = 1; # disable stdout buffering

# get latest entries.idx
my $PDB_ENTRIES = "entries.idx";
unlink $PDB_ENTRIES;
system("wget ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/$PDB_ENTRIES");

# get latest pdb_revdat.txt if it doesn't exist
my $pdbrevdat = "$Bin/pdb_revdat.txt";
if (!-s $pdbrevdat) {
	system("wget http://robetta.bakerlab.org/downloads/databases/pdb_revdat.txt.gz");
	system("gunzip pdb_revdat.txt.gz");
	system("mv pdb_revdat.txt $pdbrevdat") if (!-s $pdbrevdat);
}

# http location of rcsb pdb files
my $RCSB_PDB = "http://www.rcsb.org/pdb/files";

# optional path to hierarchical PDB directory from RCSB
my $LAB_PDB = $ENV{'PDB_DIR'} if (exists $ENV{'PDB_DIR'});
#system("rsync -auvz --port=33444 rsync.wwpdb.org::ftp/data/structures/divided/pdb/ $LAB_PDB > rsync_divided_pdb.log 2>/dev/null");

# current date
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900;
$mon++;
my $current_date_str = $mon."-".$mday."-".$year;

my (%pdbrevdat, %entries);
# get current revdat pdb codes
open(F, "$pdbrevdat") or die "ERROR! cannot open $pdbrevdat: $!\n";
while (<F>) {
	if (/^(\S{4})\s+/) {
		$pdbrevdat{$1}++;
	}
}
close(F);

# get PDB entries
open(F, $PDB_ENTRIES) or die "ERROR! cannot open $PDB_ENTRIES: $!\n";
while (<F>) {
  my @cols = split(/\t/, $_);
  if ($cols[2] =~ /^(\d+)\/(\d+)\/(\d+)$/) {
    my $pdbid = lc $cols[0];
		next if ($pdbrevdat{$pdbid});
    $entries{$pdbid} = $cols[2];
  }
}
close(F);

print "cp $pdbrevdat $pdbrevdat.$current_date_str\n";
system("cp $pdbrevdat $pdbrevdat.$current_date_str");
print "cp $pdbrevdat $pdbrevdat.bak\n";
system("cp $pdbrevdat $pdbrevdat.bak"); # backup

open(NF, ">>$pdbrevdat.$current_date_str") or die "ERROR! cannot open $pdbrevdat.$current_date_str: $!\n";
foreach my $code (keys %entries) {
	next if ($code !~ /^\S{4}$/);
	my @revdat;
	my $dir = substr($code, 1, 2);
	if (-s "$LAB_PDB/$dir/$code.pdb") {
		@revdat = &get_REVDAT( $code, "$LAB_PDB/$dir/$code.pdb" );
	} else {
		system("wget $RCSB_PDB/$code.pdb");
		if (!-s "$code.pdb") {
			warn "WARNING! cannot get REVDAT since $code.pdb does not exist at $RCSB_PDB/$code.pdb.\n";
			next;
		} else {
			@revdat = &get_REVDAT( $code, "$code.pdb" );
		}
	}
	if (!scalar@revdat) {
		warn "WARNING! REVDAT missing in $code\n";
	} else {
		print "New REVDAT entry: $code\n";
		print NF @revdat;
	}
	unlink "$code.pdb";
}
close(NF);
print "mv $pdbrevdat.$current_date_str $pdbrevdat\n";
system("mv $pdbrevdat.$current_date_str $pdbrevdat");




sub get_REVDAT {
	my ($code, $pdb) = @_;
	my @revdat;
	open(F, $pdb) or die "ERROR! cannot open $pdb: $!\n";
	while (<F>) {
		# get release date only
		if (/^REVDAT\s+\d+\s+\d+\-\w\w\w\-\d+\s+\S{4}\s+(\d+)/ && $1 eq '0') {
			push(@revdat, "$code $_");
		}
	}
	close(F);
	return @revdat;
}

