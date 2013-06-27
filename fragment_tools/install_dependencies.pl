#!/usr/bin/perl -w
use strict;
use FindBin qw( $Bin );
$| = 1; # disable stdout buffering

if (!scalar@ARGV || $ARGV[0] !~ /^(standard|overwrite)$/) {
	print "\n";
	print "USAGE: $0 <standard|overwrite> [skip_nr]\n\n";
	print "This script installs blast, psipred, sparks-x, and the NCBI non-redundant (nr) database.\n";
	print "The <standard> option will only install what is missing.\n";
	print "The <overwrite> option will do a fresh installation.\n\n";
	print "Install locations: \n";
	print " $Bin/blast\n";
	print " $Bin/psipred\n";
	print " $Bin/sparks-x\n";
	print " $Bin/databases/nr\n";
	print " $Bin/databases/nr_pfilt\n";
	print "\n";
	print "Requires ~73+ Gigs of free disk space.\n";
	print "The NCBI sequence databases are very large.\n";
	print "Please be patient when running this script.\n";
	print "\n";
	print "FOR ACADEMIC USE ONLY.\n\n";
	exit(1);
}
my $overwrite = ($ARGV[0] eq 'overwrite') ? 1 : 0;
my @packages_to_clean;

my $skip_nr = 0;
foreach my $arg (@ARGV) {
	$skip_nr = 1 if ($arg =~ /^skip_nr\s*$/);
}

chdir($Bin);

# blast binaries
# from ftp://ftp.ncbi.nih.gov/blast/executables/release/2.2.17/
if ($overwrite || !-d "$Bin/blast/bin" || !-d "$Bin/blast/data") {
	# try to figure out what package to install
	my $package = "blast-2.2.17-ia32-linux.tar.gz";
	my $proc = `uname -m`;
	if ($proc =~ /x86_64/) {
		$package = "blast-2.2.17-x64-linux.tar.gz";
	} elsif ($proc =~ /ia64/) {
		$package = "blast-2.2.17-ia64-linux.tar.gz";
	}
	my $url = "ftp://ftp.ncbi.nih.gov/blast/executables/release/2.2.17/$package";
	print "INSTALLING BLAST from $url ....\n";
	system("rm -rf blast") if (-d "blast");  # clean up interrupted attempts
	system("wget -N $url");
	system("tar -zxvf $package");
	push(@packages_to_clean, "$Bin/$package");
	system("mv blast-2.2.17 blast");
	(-d "$Bin/blast/bin" && -d "$Bin/blast/data") or die "ERROR! blast installation failed!\n";
}
my $blast_credit =<<BLASTCREDIT;

>>>> PLEASE CITE THE FOLLOWING FOR BLAST <<<<

Altschul, S.F., Madden, T.L., Schäffer, A.A., Zhang, J., Zhang, Z., Miller
W. & Lipman, D.J (1997) "Gapped BLAST and PSI-BLAST: a new generation of
protein database search programs." Nucleic Acids Res. 25:3389-3402.

BLASTCREDIT
chdir($Bin);

# psipred
# from http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/psipred3.3.tar.gz
if ($overwrite || !-d "$Bin/psipred/bin" || !-d "$Bin/psipred/data") {
	my $package = "psipred3.3.tar.gz";
	my $url = "http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/$package";
	print "INSTALLING PSIPRED from $url ....\n";
	system("rm -rf $Bin/psipred") if (-d "$Bin/psipred"); # clean up interrupted attempts
	(mkdir("$Bin/psipred")) or die "ERROR! cannot mkdir $Bin/psipred: $!\n";
	chdir("$Bin/psipred/");
	system("wget -N $url");
	if (!-s $package ) {
		$url = "http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/old/$package";
		print "wget failed, trying $url ....\n";
		system("wget -N $url");
	}
	system("tar -zxvf $package");
	push(@packages_to_clean, "$Bin/psipred/$package");
	chdir("$Bin/psipred/src");
	system("make"); # build from src
	my @binaries = qw( chkparse pfilt psipass2 psipred seq2mtx );
	foreach my $binary (@binaries) {
		system("mv $binary $Bin/psipred/bin/");
	}
	(-d "$Bin/psipred/bin" && -d "$Bin/psipred/data") or die "ERROR! psipred installation failed!\n";
}
my $psipred_credit =<<PSIPREDCREDIT;

>>>> PLEASE CITE THE FOLLOWING FOR PSIPRED <<<<

Jones, D.T. (1999) Protein secondary structure prediction based on
position-specific scoring matrices. J. Mol. Biol. 292:195-202.

PSIPREDCREDIT
chdir($Bin);

# SPARKS-X/SPINE-X
# from http://sparks.informatics.iupui.edu/yueyang/download/SPARKS-X/sparksx-1.tgz
if ($overwrite || !-d "$Bin/sparks-x/bin" || !-d "$Bin/sparks-x/data") {
	my $package = "sparksx-1.tgz";
	my $url = "http://sparks.informatics.iupui.edu/yueyang/download/SPARKS-X/$package";
	print "INSTALLING SPARKS-X from $url ....\n";
	system("rm -rf sparks-x") if (-d "sparks-x"); # clean up interrupted attempts
	system("wget -N $url");
	system("tar -zxvf $package");
	push(@packages_to_clean, "$Bin/$package");
	# update paths to sparks-x directory in sparks-x scripts
	# instead of using the SPARKSXDIR environment variable.
	my @scripts = qw( SPINE-X/spineX.pl bin/buildinp_query.sh bin/psiblast.sh bin/scan1.sh bin/scan_multi.sh );
	foreach my $script (@scripts) {
		system("cp $Bin/sparks-x/$script $Bin/sparks-x/$script.orig");
		open(F, "$Bin/sparks-x/$script.orig") or die "ERROR! cannot open $Bin/sparks-x/$script.orig: $!\n";
		open(NF, ">$Bin/sparks-x/$script") or die "ERROR! cannot open $Bin/sparks-x/$script: $!\n";
		while (my $l = <F>) {
			if ($l =~ /^spxdir=/) {
				print NF "#".$l;
				print NF "# the following 2 lines were added by Rosetta/tools/fragment_tools/install_dependencies.pl\n";
				print NF "spxdir=\$(dirname \$(readlink -f \$0))\n";
				print NF "spxdir=\${spxdir\%/*}\n";
				next;
			}
			if ($l !~ /^\s*#/ && $l =~ /[\`"']\s*cp\s+\-n\s+/) { # replace cp -n with cp since -n option may not exist
				my $nl = $l;
				$nl =~ s/([\`"']\s*)cp\s+\-n\s+/$1cp /gs;
				chomp $nl;
				print NF "#".$l;
				print NF "# the following line was added by Rosetta/tools/fragment_tools/install_dependencies.pl\n";
				print NF "$nl\n";
				next;
			}
			print NF $l;
		}
		close(F);
		close(NF);
	}
	chdir("$Bin/sparks-x");
	system("ln -sf ../blast ./");
	system("ln -sf ../databases/ ./blast-NR");
	(-d "$Bin/sparks-x/bin" && -d "$Bin/sparks-x/data") or die "ERROR! sparks-x installation failed!\n";
}
my $sparksx_credit =<<SPARKSXCREDIT;

>>>> PLEASE CITE THE FOLLOWING FOR SPARKS-X <<<<

Y. Yang, E. Faraggi, H. Zhao and Y. Zhou. Improving protein fold recognition
and template-based modeling by employing probabilistic-based matching between
predicted one-dimensional structural properties of the query and corresponding
native properties of templates. Bioinformatics 27, 2076-2082 (2011)

SPARKSXCREDIT
chdir($Bin);

# clean up
foreach my $pkg (@packages_to_clean) { unlink $pkg; }

our $datdir = "$Bin/databases";
(-d $datdir || mkdir($datdir)) or die "ERROR! cannot mkdir $datdir: $!\n";

# may want to add options here for lighter weight sequence databases like:
#  ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz  90% clustering
#  ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz  50% clustering
#
#  Sequence searches would be much faster

# NR database
if (!$skip_nr && ($overwrite || !-s "$datdir/nr.pal")) {
	chdir($datdir);
	system("wget -N http://www.ncbi.nlm.nih.gov/blast/docs/update_blastdb.pl");
	die "ERROR! wget http://www.ncbi.nlm.nih.gov/blast/docs/update_blastdb.pl failed.\n" if (!-s "update_blastdb.pl");
	print "Fetching NR database from NCBI. Be very patient ......\n";
	system("rm nr.*"); # clean up interrupted attempts
	$SIG{INT} = \&clean_nr_tgz;
	(system("perl update_blastdb.pl nr") == 0) or do { &clean_nr_tgz; };
	$SIG{INT} = \&clean_nr;
	foreach my $f (glob("$datdir/nr.*tar.gz")) {
		my $unzipped = $f;
		$unzipped =~ s/\.tar\.gz$//;
		if (!-s $unzipped) {
			(system("tar -zxvf $f") == 0) or do { &clean_nr; };
			unlink $f; # save disk space
		}
	}
	$SIG{INT} = 'DEFAULT';
	(-s "$datdir/nr.pal") or die "ERROR! nr database installation failed!\n";
}
chdir($Bin);

# p_filt NR database
if (!$skip_nr && ($overwrite || !-s "$datdir/nr_pfilt.pal")) {
	chdir($datdir);
	print "Generating nr fasta. Be very very patient ......\n";
	my $cmd = "$Bin/blast/bin/fastacmd -D 1 > nr";
	print "$cmd\n";
	(system($cmd) == 0) or die "ERROR! $cmd failed.\n";
	print "Generating nr_pfilt fasta. Be very very very patient ......\n";
	$cmd = "$Bin/psipred/bin/pfilt nr > nr_pfilt";
	print "$cmd\n";
	(system($cmd) == 0) or die "ERROR! $cmd failed.\n";
	print "Formatting nr_pfilt fasta. Be very very very very patient ......\n";
	system("rm nr_pfilt.*"); # clean up interrupted attempts
	$SIG{INT} = \&clean_nr_pfilt;
	(system("$Bin/blast/bin/formatdb -o T -i nr_pfilt") == 0) or do { &clean_nr_pfilt; };
	$SIG{INT} = 'DEFAULT';
	(-s "$datdir/nr_pfilt.pal") or die "ERROR! nr_pfilt database installation failed!\n";
}

print "\n";
print " Installed\n";
print "     blast: $Bin/blast\n";
print "   psipred: $Bin/psipred\n";
print "  sparks-x: $Bin/sparks-x\n";
print "        nr: $Bin/databases/nr\n" if (!$skip_nr);
print "  nr_pfilt: $Bin/databases/nr_pfilt\n" if (!$skip_nr);
print "\nDone!\n";

print $blast_credit;
print $psipred_credit;
print $sparksx_credit;

print "FOR ACADEMIC USE ONLY.\n\n";

## catch ctrl-c
sub clean_nr_tgz {
	$SIG{INT} = \&clean_nr_tgz;
	warn "Aborting....\n";
	system("rm $datdir/nr.*tar.gz");
	die "Aborted!\n";
}
sub clean_nr {
	$SIG{INT} = \&clean_nr;
	warn "Aborting....\n";
	system("rm $datdir/nr.??.p*");
	die "Aborted!\n";
}
sub clean_nr_pfilt {
	$SIG{INT} = \&clean_nr_pfilt;
	warn "Aborting....\n";
	system("rm $datdir/nr_pfilt.??.p*");
	die "Aborted!\n";
}


