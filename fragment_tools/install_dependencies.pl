#!/usr/bin/perl

my $pwd=`pwd`;
chomp $pwd;

print "\n";
print "This script installs blast, psipred, sparks-x, and the ncbi non-redundant (nr) database.\n";
print "Requires ~73+ Gigs of free disk space.\n"; 

# blast binaries
# from ftp://ftp.ncbi.nih.gov/blast/executables/release/2.2.17/
if (!-d "$pwd/blast/bin" || !-d "$pwd/blast/data") {
	my $package = "blast-2.2.17-ia32-linux.tar.gz";
	my $proc = `uname -m`;
	if ($proc =~ /x86_64/) {
		$package = "blast-2.2.17-x64-linux.tar.gz";
	} elsif ($proc =~ /ia64/) {
		$package = "blast-2.2.17-ia64-linux.tar.gz";
	}
	my $url = "ftp://ftp.ncbi.nih.gov/blast/executables/release/2.2.17/$package";
	print "INSTALLING BLAST from $url ....\n";
	system("wget $url");
	system("tar -zxvf $package");
	unlink "$package";
	system("mv blast-2.2.17 blast");
	(-d "$pwd/blast/bin" && -d "$pwd/blast/data") or die "ERROR! blast installation failed!\n";
}
chdir($pwd);

# psipred
# from http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/psipred3.3.tar.gz
if (!-d "$pwd/psipred/bin" || !-d "$pwd/psipred/data") {
	my $package = "psipred3.3.tar.gz";
	my $url = "http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/$package";
	print "INSTALLING PSIPRED from $url ....\n";
	system("wget $url");
	if (!-s $package ) {
		$url = "http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/old/$package";
 		print "wget failed, trying $url ....\n";
		system("wget $url");
	}
	(-d "$pwd/psipred" || mkdir("$pwd/psipred")) or die "ERROR! cannot mkdir $pwd/psipred: $!\n";
	system("mv $pwd/$package $pwd/psipred/");
	chdir("$pwd/psipred/");
	system("tar -zxvf $package");
	unlink $package;
	chdir("$pwd/psipred/src");
	system("make"); # build from src
	my @binaries = qw( chkparse pfilt psipass2 psipred seq2mtx );
	foreach my $binary (@binaries) {
		system("mv $binary $pwd/psipred/bin/");
	}
	(-d "$pwd/psipred/bin" && -d "$pwd/psipred/data") or die "ERROR! psipred installation failed!\n";
}
my $psipred_credit =<<PSIPREDCREDIT;

>>>> PLEASE CITE THE FOLLOWING FOR PSIPRED <<<<

Jones, D.T. (1999) Protein secondary structure prediction based on
position-specific scoring matrices. J. Mol. Biol. 292:195-202.
PSIPREDCREDIT
print $psipred_credit;
chdir($pwd);

# SPARKS-X/SPINE-X
# from http://sparks.informatics.iupui.edu/yueyang/download/SPARKS-X/sparksx-1.tgz
if (!-d "$pwd/sparks-x/bin" || !-d "$pwd/sparks-x/data") {
	my $package = "sparksx-1.tgz";
	my $url = "http://sparks.informatics.iupui.edu/yueyang/download/SPARKS-X/$package";
	print "INSTALLING SPARKS-X from $url ....\n";
	system("wget $url");
	system("tar -zxvf $package");
	unlink $package;
	my @scripts = qw( buildinp_query.sh psiblast.sh scan1.sh scan_multi.sh );
	foreach my $script (@scripts) {
		system("cp $pwd/sparks-x/bin/$script $pwd/sparks-x/bin/$script.orig");
		open(F, "$pwd/sparks-x/bin/$script.orig") or die "ERROR! cannot open $pwd/sparks-x/bin/$script.orig: $!\n";
		open(NF, ">$pwd/sparks-x/bin/$script") or die "ERROR! cannot open $pwd/sparks-x/bin/$script: $!\n";
		while (my $l = <F>) {
			if ($l =~ /^spxdir=/) {
				print NF "spxdir=$pwd/sparks-x\n";
				next;
			} elsif ($l =~ /buildinp\.py/) {
				next; # skip this step
			}
			print NF $l;
		}
		close(F);
		close(NF);
	}
	chdir("$pwd/sparks-x");
	system("ln -sf ../blast ./");
	system("ln -sf ../databases/ ./blast-NR");
	(-d "$pwd/sparks-x/bin" && -d "$pwd/sparks-x/data") or die "ERROR! sparks-x installation failed!\n";
}
my $sparksx_credit =<<SPARKSXCREDIT;

>>>> PLEASE CITE THE FOLLOWING FOR SPARKS-X <<<<

Y. Yang, E. Faraggi, H. Zhao and Y. Zhou. Improving protein fold recognition
and template-based modeling by employing probabilistic-based matching between
predicted one-dimensional structural properties of the query and corresponding
native properties of templates. Bioinformatics 27, 2076-2082 (2011)

SPARKSXCREDIT
print $sparksx_credit;
chdir($pwd);

# NR database
my $datdir = "$pwd/databases";
(-d $datdir || mkdir($datdir)) or die "ERROR! cannot mkdir $datdir: $!\n";
if (!-s "$datdir/nr.pal") {
	chdir($datdir);
	if (!-s "update_blastdb.pl") {
		system("wget http://www.ncbi.nlm.nih.gov/blast/docs/update_blastdb.pl");
	}
	print "Fetching NR database from NCBI. Be very patient ......\n";
	$SIG{INT} = \&clean_nr_tgz; # make sure we cleanup partially downloaded files if CTRL+C'd
	system("perl update_blastdb.pl nr");
	$SIG{INT} = \&clean_nr; # make sure we cleanup partially downloaded files if CTRL+C'd
	foreach my $f (glob("$datdir/nr.*tar.gz")) {
		my $unzipped = $f;
		$unzipped =~ s/\.tar\.gz$//;
		system("tar -zxvf $f") if (!-s $unzipped);
	}
	(-s "$datdir/nr.pal") or die "ERROR! nr database installation failed!\n";
	system("rm $datdir/nr.*tar.gz"); # save disk space
}
chdir($pwd);

# p_filt NR database
if (!-s "$datdir/nr_pfilt.pal") {
	chdir($datdir);
	if (!-s "nr") {
		print "Generating nr fasta. Be very very patient ......\n";
		$SIG{INT} = \&clean_nr_fasta;
		print "$pwd/blast/bin/fastacmd -D 1 > nr\n";
		system("$pwd/blast/bin/fastacmd -D 1 > nr");
	}
	if (!-s "nr_pfilt") {
		print "Generating nr_pfilt fasta. Be very very very patient ......\n";
		$SIG{INT} = \&clean_nr_pfilt_fasta;
		print "$pwd/psipred/bin/pfilt nr > nr_pfilt\n";
		system("$pwd/psipred/bin/pfilt nr > nr_pfilt");
	}
	$SIG{INT} = \&clean_nr_pfilt;
	print "Formatting nr_pfilt fasta. Be very very very very patient ......\n";
	system("$pwd/blast/bin/formatdb -o T -i nr_pfilt");
	(-s "$datdir/nr_pfilt.pal") or die "ERROR! nr_pfilt database installation failed!\n";
}

print " Installed\n";
print "     blast: $pwd/blast\n";
print "   psipred: $pwd/psipred\n";
print "  sparks-x: $pwd/sparks-x\n";
print "        nr: $pwd/databases/nr\n";
print "  nr_pfilt: $pwd/databases/nr_pfilt\n";
print "\nDone!\n";


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
sub clean_nr_fasta {
	$SIG{INT} = \&clean_nr_fasta;
	warn "Aborting....\n";
	system("rm $datdir/nr");
	die "Aborted!\n";
}
sub clean_nr_pfilt_fasta {
  $SIG{INT} = \&clean_nr_pfilt_fasta;
  warn "Aborting....\n";
  system("rm $datdir/nr_pfilt");
  die "Aborted!\n";
}
sub clean_nr_pfilt {
	$SIG{INT} = \&clean_nr_pfilt;
	warn "Aborting....\n";
	system("rm $datdir/nr_pfilt.??.p*");
	die "Aborted!\n";
}


