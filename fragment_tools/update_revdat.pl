#!/usr/bin/perl
use FindBin qw( $Bin );
$| = 1; # disable stdout buffering

# This script grabs the latest Robetta specific file 'pdb_revdat.txt.gz'
# from the Robetta server, which updates the file weekly after the
# RCSB updates the PDB on Wednesdays +0 UTC.

# pdb_revdat.txt is used in make_fragments.pl to date filter for benchmarking

my $INTERNET_HOST = "";
if (scalar@ARGV) {
	$INTERNET_HOST = shift@ARGV;
	chomp $INTERNET_HOST;
}

# optional path to hierarchical PDB directory from RCSB
my $LAB_PDB = $ENV{'PDB_DIR'} if (exists $ENV{'PDB_DIR'});
#system("rsync -auvz --port=33444 rsync.wwpdb.org::ftp/data/structures/divided/pdb/ $LAB_PDB > rsync_divided_pdb.log 2>/dev/null");

# http location of rcsb pdb files
my $RCSB_PDB = "http://www.rcsb.org/pdb/files";

# current date
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900;
$mon++;
my $current_date_str = $year.sprintf("%2.2d", $mon).sprintf("%2.2d",$mday);

# get latest pdb_revdat.txt if it doesn't exist
my $pdbrevdat = "$Bin/pdb_revdat.txt";
if (!-s $pdbrevdat || (scalar@ARGV && $ARGV[0] eq 'overwrite')) {
	if ($INTERNET_HOST) {
		system("ssh $INTERNET_HOST 'cd $Bin; wget -N http://robetta.bakerlab.org/downloads/databases/pdb_revdat.txt.gz'");
	} else {
		system("cd $Bin; wget -N http://robetta.bakerlab.org/downloads/databases/pdb_revdat.txt.gz");
	}
	system("gunzip -c $Bin/pdb_revdat.txt.gz > $pdbrevdat.$current_date_str");
	system("mv $pdbrevdat.$current_date_str $pdbrevdat");
} else {
	my $epoch_timestamp = (stat("$Bin/pdb_revdat.txt.gz"))[9];
	# Subtract an hour to account for time Robetta may take to update the file on Wednesdays.
	# There's a time window when the RCSB will be updated before pdb_revdat.txt.gz, but this
	# is most likely never going to be an issue when date filtering for benchmarking.
	my $check_timestamp = $epoch_timestamp-3600;
	my $current_epoch = time;
	my $diff =  ($current_epoch-$check_timestamp)/604800;
	if ($diff > 1) {
		# older than a week so lets try to update
		if ($INTERNET_HOST) {
			system("ssh $INTERNET_HOST 'cd $Bin; wget -N http://robetta.bakerlab.org/downloads/databases/pdb_revdat.txt.gz'");
		} else {
			system("cd $Bin; wget -N http://robetta.bakerlab.org/downloads/databases/pdb_revdat.txt.gz");
		}
		my $new_epoch_timestamp = (stat("$Bin/pdb_revdat.txt.gz"))[9];
		if ($new_epoch_timestamp > $epoch_timestamp) {
			system("gunzip -c $Bin/pdb_revdat.txt.gz > $pdbrevdat.$current_date_str");
			system("mv $pdbrevdat.$current_date_str $pdbrevdat");
		}
	}
}

if (!scalar@ARGV || $ARGV[0] ne 'force') {
	exit(0);
}

# If called with <force>, download new PDBs from the RCSB if they exist and add them to pdb_revdat.txt

# get latest entries.idx
my $PDB_ENTRIES = "entries.idx";
unlink $PDB_ENTRIES;
system("wget ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/$PDB_ENTRIES");

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

