#!/usr/bin/perl

use FindBin qw($Bin);
use Fcntl ':flock';

# The PDB archive is updated on Wednesday 00:00 UTC which is Tuesday 4:00 PM PST
# Return if it's not Tuesday-Thursday (to cover all timezones)
my $log = "$Bin/logs";
if (!scalar@ARGV && -s $log) {
  my $wday = (localtime(time))[6]; # 0 is Sunday, 1 is Monday, 2 is Tuesday ....
  exit(0) if ($wday < 2 || $wday > 4);
}

# run rsyncPDB.sh with a minimum interval so we don't rsync the RCSB too often
my $min_interval = 60 * 30; # 30 minutes

my $ss_dis = "http://www.rcsb.org/pdb/files/ss_dis.txt.gz";
my $lockfile = "$Bin/update_rcsb_data.lock";

# get an exclusive lock
open(LOCKF, ">>$lockfile") or die "Could not open '$lockfile' - $!";
flock(LOCKF, LOCK_EX) or die "Could not lock '$lockfile' - $!";

my $epoch = time();

## get the age of the rsyncPDB.sh log file
my $diff = (-s $log) ? $epoch - (stat($log))[8] :  $min_interval+1;

## get the size of the last rsyncPDB.sh update
my $prevsize = &gettotalsizefromlog($log);
# rsync PDB data if within the minimum interval
# rsyncPDB.sh will only update if the RCSB archive was updated
if ($diff > $min_interval) {
  print "$Bin/rsyncPDB.sh $Bin\n";
  (system("$Bin/rsyncPDB.sh $Bin") == 0) or warn "WARNING! $Bin/rsyncPDB.sh $Bin failed\n";
  ## update timestamp of rsyncPDB.sh log file
  (system("touch $log") == 0) or &handle_error("ERROR! touch $log failed");
}
# get the size of the recent update
my $currsize = &gettotalsizefromlog($log);

# only update ss_dis.txt if the RCSB archive was updated and ss_dis.txt is older than 4 days
# just to be safe so we don't update ss_dis.txt too often. ss_dis.txt is a big file and can
# only be retrieved via http and the RCSB http header does not provide Last-Modified info.
my $ss_dis_diff = (-s "$Bin/ss_dis.txt") ? $epoch - (stat("$Bin/ss_dis.txt"))[8] : 7 * 60 * 60 * 24;
if (  $prevsize != $currsize ||            # update if rsyncPDB.sh updated files
      $ss_dis_diff > 4 * 60 * 60 * 24  ) { # update if older than 4 days
  print "UPDATING ss_dis.txt last prev rsync size: $prevsize curr rsync size: $currsize\n";
  my $tmpdir = "$Bin/tmp";
  (-d $tmpdir || mkdir($tmpdir)) or &handle_error("ERROR! cannot mkdir $tmpdir: $!");
  # get ss_dis.txt from RCSB
  (system("wget $ss_dis -O $tmpdir/ss_dis.txt.gz; gunzip -f $tmpdir/ss_dis.txt.gz") == 0) or
    &handle_error("ERROR! wget $ss_dis failed");
  (system("mv $tmpdir/ss_dis.txt $Bin/") == 0) or &handle_error("ERROR! mv $tmpdir/ss_dis.txt $Bin/ failed");
  (system("touch $Bin/ss_dis.txt") == 0) or &handle_error("ERROR! touch $Bin/ss_dis.txt failed");
}
close(LOCKF) or die "Could not write '$lockfile' - $!";


sub handle_error {
  my $msg = shift;
	warn "$msg\n";
	close(LOCKF) or die "Could not close '$lockfile' - $!";
	exit(1);
}

sub gettotalsizefromlog {
  my $l = shift;
  my $size = 0;
  return 0 if (!-f $l);
  open(F, $l) or die "ERROR! cannot open $l: $!\n";
  while (<F>) {
	  if (/total\s+size\s+is\s+(\d+)/) {
		  $size = $1;
    }
  }
	close(F);
	return $size;
}


