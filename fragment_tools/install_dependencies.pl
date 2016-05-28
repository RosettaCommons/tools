#!/usr/bin/perl -w
use strict;
use FindBin qw( $Bin );
$| = 1; # disable stdout buffering

my $installtype = "";
foreach my $arg (@ARGV) {
	if ($arg =~ /^(standard|overwrite)\s*$/) {
		$installtype = $1;
	}
	if ($arg !~ /^(standard|overwrite|nr|uniref90|uniref50|skip_nr)\s*$/) {
		$installtype = "";
		last;
	}
}
if (!$installtype) {
	print "\n";
	print "USAGE: $0 <standard|overwrite> [nr(default)|uniref90|uniref50|skip_nr]\n\n";
	print "This script installs blast, psipred, csblast, sparks-x, and the NCBI non-redundant (nr) database.\n";
	print "<standard> will only install what is missing.\n";
	print "<overwrite> will do a fresh installation.\n";
	print "[uniref90] will use the uniref90 sequence database (NCBI nr is the default.)\n";
	print "[uniref50] will use the uniref50 sequence database.\n";
	print "[skip_nr] skips the sequence database installation.\n";
	print "uniref databases will OVERWRITE previously installed nr if it exists\n\n";
	print "Install locations: \n";
	print " $Bin/blast\n";
	print " $Bin/psipred\n";
	print " $Bin/csblast\n";
	print " $Bin/sparks-x\n";
	print " $Bin/databases/nr\n";
	print " $Bin/databases/nr_pfilt\n";
	print "\n";
	print "Requires ~73+ Gigs of free disk space for the NCBI nr database.\n";
	print "The sequence database is very large.\n";
	print "Please be patient when running this script.\n";
	print "\n";
	print "FOR ACADEMIC USE ONLY.\n\n";
	exit(1);
}
my $overwrite = ($installtype eq 'overwrite') ? 1 : 0;
my @packages_to_clean;

my $skip_nr = 0;
my $database = "nr";
foreach my $arg (@ARGV) {
	$skip_nr = 1 if ($arg =~ /^skip_nr\s*$/);
	if ($arg =~ /^(uniref90|uniref50)\s*$/) { $database = $1; }
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
	my $url = "ftp://ftp.ncbi.nih.gov/blast/executables/legacy/2.2.17/$package";
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
		$url = "http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/old_versions/$package";
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

# csblast https://github.com/cangermueller/csblast
if ($overwrite || !-d "$Bin/csblast/bin") {
	my $package = "v2.2.3.tar.gz";
	my $url = "https://github.com/cangermueller/csblast/archive/$package";

	print "INSTALLING CSBLAST from $url\n";

	system("rm -rf $Bin/csblast") if (-d "$Bin/csblast"); # clean up interrupted attempts
	(mkdir("$Bin/csblast")) or die "ERROR! cannot mkdir $Bin/csblast$!\n";
	chdir("$Bin/csblast/");
	system("wget -N $url");
	system("tar -zxvf $package");
  my $sph_package = "sparsehash-2.0.3.tar.gz";
  system("wget -N https://github.com/sparsehash/sparsehash/archive/$sph_package");
	system("tar -zxvf $sph_package");
	push(@packages_to_clean, "$Bin/csblast/$package", "$Bin/csblast/$sph_package");

	chdir("$Bin/csblast/sparsehash-sparsehash-2.0.3");
  system("./configure --prefix=$Bin/csblast/local");
	system("make install"); # build from src

	chdir("$Bin/csblast");
  system("mv csblast-*/* .");
	chdir("$Bin/csblast/src");
  system('sed -i -e "s|^INC.*|INC = -I../local/include|" Makefile');
  system("make csblast csbuild");

	(-d "$Bin/csblast/bin") or die "ERROR! psipred installation failed!\n";
}

my $csblast_credit =<<CSBLASTCREDIT;

>>>> PLEASE CITE THE FOLLOWING FOR CSBLAST <<<<

Angermüller, C., Biegert, A., & Söding, J. (2012). Discriminative modelling of
context-specific amino acid substitution probabilities. Bioinformatics, 28(24),
3240-3247.

Biegert, A., & Söding, J. (2009). Sequence context-specific profiles for
homology searching.  Proceedings of the National Academy of Sciences,106(10),
3770-3775.

CSBLASTCREDIT
chdir($Bin);

# SPARKS-X/SPINE-X
# from http://sparks-lab.org/pmwiki/download/yueyang/SPARKS-X/sparksx-1.tgz
if ($overwrite || !-d "$Bin/sparks-x/bin" || !-d "$Bin/sparks-x/data") {
	my $package = "sparksx-1.tgz";
	my $url = "http://sparks-lab.org/pmwiki/download/yueyang/SPARKS-X/$package";
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

# options for lighter weight sequence databases
#  ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz  90% clustering
if (!$skip_nr && $database eq "uniref90" && ($overwrite || !-f "$datdir/uniref90")) {
	chdir($datdir);
	print "Fetching uniref90 database. Be very patient ......\n";
	system("rm $datdir/uniref90* $datdir/nr*"); # clean up interrupted attempts
	$SIG{INT} = \&clean_uniref90;
	(system("wget -N ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz") == 0) or
		do { &clean_uniref90; };
	print "Extracting uniref90.fasta.gz to nr fasta. Be very very patient ......\n";
	(system("gunzip $datdir/uniref90.fasta.gz") == 0) or do { &clean_uniref90; };
	(system("mv $datdir/uniref90.fasta $datdir/nr") == 0) or do { &clean_uniref90; };
	print "Formating nr fasta. Be very very very patient ......\n";
	(system("$Bin/blast/bin/formatdb -o T -i $datdir/nr") == 0) or do { &clean_uniref90; };
	$SIG{INT} = 'DEFAULT';
	(-s "$datdir/nr.pal") or die "ERROR! uniref90 as $datdir/nr database installation failed!\n";
	(system("ln -s nr uniref90") == 0) or do { &clean_uniref90; }; # create a link just to know the nr is uniref
}
chdir($Bin);

#  ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz  50% clustering
if (!$skip_nr && $database eq "uniref50" && ($overwrite || !-s "$datdir/uniref50")) {
	chdir($datdir);
	print "Fetching uniref50 database. Be very patient ......\n";
	system("rm $datdir/uniref50* $datdir/nr*"); # clean up interrupted attempts
	$SIG{INT} = \&clean_uniref50;
	(system("wget -N ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz") == 0) or
		do { &clean_uniref50; };
	print "Extracting uniref50.fasta.gz to nr fasta. Be very very patient ......\n";
	(system("gunzip $datdir/uniref50.fasta.gz") == 0) or do { &clean_uniref50; };
	(system("mv $datdir/uniref50.fasta $datdir/nr") == 0) or do { &clean_uniref50; };
	print "Formating nr fasta. Be very very very patient ......\n";
	(system("$Bin/blast/bin/formatdb -o T -i $datdir/nr") == 0) or do { &clean_uniref50; };
	$SIG{INT} = 'DEFAULT';
	(-s "$datdir/nr.pal") or die "ERROR! uniref50 as $datdir/nr database installation failed!\n";
	(system("ln -s nr uniref50") == 0) or do { &clean_uniref50; }; # create a link just to know the nr is uniref
}
chdir($Bin);

# big fat NR database!
if (!$skip_nr && $database eq "nr" && ($overwrite || !-s "$datdir/nr.pal")) {
	chdir($datdir);
	print "Fetching NR database from NCBI. Be very patient ......\n";
	system("wget -N http://www.ncbi.nlm.nih.gov/blast/docs/update_blastdb.pl");
	die "ERROR! wget http://www.ncbi.nlm.nih.gov/blast/docs/update_blastdb.pl failed.\n" if (!-s "$datdir/update_blastdb.pl");
	system("rm $datdir/nr*"); # clean up interrupted attempts
	$SIG{INT} = \&clean_nr_tgz;
	(system("perl $datdir/update_blastdb.pl nr") == 0) or do { &clean_nr_tgz; };
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
	(-s "$datdir/nr.pal") or die "ERROR! $datdir/nr database installation failed!\n";
}
chdir($Bin);

# p_filt NR database from installed nr
if (!$skip_nr && ($overwrite || !-s "$datdir/nr_pfilt.pal")) {
	chdir($datdir);
	print "Generating nr fasta. Be very patient ......\n";
	my $cmd = "$Bin/blast/bin/fastacmd -D 1 > $datdir/nr";
	(system($cmd) == 0) or die "ERROR! $cmd failed.\n";
	print "Generating nr_pfilt fasta. Be very very patient ......\n";
	$cmd = "$Bin/psipred/bin/pfilt $datdir/nr > $datdir/nr_pfilt";
	(system($cmd) == 0) or die "ERROR! $cmd failed.\n";
	print "Formatting nr_pfilt fasta. Be very very very patient ......\n";
	system("rm $datdir/nr_pfilt.*"); # clean up interrupted attempts
	$SIG{INT} = \&clean_nr_pfilt;
	(system("$Bin/blast/bin/formatdb -o T -i $datdir/nr_pfilt") == 0) or do { &clean_nr_pfilt; };
	$SIG{INT} = 'DEFAULT';
	(-s "$datdir/nr_pfilt.pal") or die "ERROR! $datdir/nr_pfilt database installation failed!\n";
}

print "\n";
print " Installed\n";
print "     blast: $Bin/blast\n";
print "   psipred: $Bin/psipred\n";
print "   csblast: $Bin/csblast\n";
print "  sparks-x: $Bin/sparks-x\n";
if (!$skip_nr) {
	print "        nr: $datdir/nr\n" if ($database eq "nr");
	print "        nr: $datdir/uniref90\n" if ($database eq "uniref90");
	print "        nr: $datdir/uniref50\n" if ($database eq "uniref50");
	print "  nr_pfilt: $datdir/nr_pfilt\n";
}
print "\nDone!\n";

print $blast_credit;
print $psipred_credit;
print $csblast_credit;
print $sparksx_credit;

print "FOR ACADEMIC USE ONLY.\n\n";

## catch ctrl-c
sub clean_uniref90 {
	$SIG{INT} = \&clean_uniref90;
	warn "Aborting....\n";
	system("rm $datdir/uniref90*");
	die "Aborted!\n";
}
sub clean_uniref50 {
	$SIG{INT} = \&clean_uniref50;
	warn "Aborting....\n";
	system("rm $datdir/uniref50*");
	die "Aborted!\n";
}
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


