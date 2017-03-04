#!/usr/bin/perl

# This script is meant to convert a rosetta++ formatted fragment file to
# a minimal fragment format which consists of a single line for each fragment
# that contains the vall line number followed by a space followed by the fragment
# length for each fragment. The very first line should start with a # followed
# by a space followed by "index" followed by the original vall starting line #,
# ending line #, previous residue key number, and name for reference. The fragment
# reader uses this line to determine which torsions only vall to use to read in
# the fragment torsions. For each blank line in the indexed fragment file, the
# position is iterated, starting at position 1 of the target sequence. The
# torsions only vall is just the vall torsions phi,psi,omega on each line.


my %opts = &getCommandLineOptions();

my $vall = $opts{vall};
my $vallname = $vall;
$vallname =~ s/^.*\/([^\/]+)\s*$/$1/gs;
$vallname =~ s/(\.gz|\.Z)\s*$//gs;

my $phicol;

my $debug = 0;

# determine torsions columns in vall
print "Reading $vall to get phi,psi,omega torsions\n";
print "  Assuming columns 1, 2, and 3 are the PDB code, chain, and seqres position.\n";

if ($vall =~ /(\.gz|\.Z)\s*$/) {
	open(F, "gzip -dc $vall |") or die "ERROR! cannot open $vall: $!\n";
} else {
	open(F, $vall) or die "ERROR! cannot open $vall: $!\n";
}
my %vall;
my $vallline = 0;
while (<F>) {
	next if (/^#/);
	my @cols = split(/\s+/,$_);
	my $tmpphicol;
	if (!defined $phicol) {
		for (my $i=0; $i<=$#cols-3; $i++) {
			next if ($cols[$i] !~ /^\-?\d{1,3}\.\d{3}$/);
			next if ($cols[$i+1] !~ /^\-?\d{1,3}\.\d{3}$/);
			next if ($cols[$i+2] !~ /^\-?\d{1,3}\.\d{3}$/);
			my $phi = eval($cols[$i]);
			my $psi = eval($cols[$i+1]);
			my $omega = eval($cols[$i+2]);
			if ($phi <= 180 && $phi >= -180 && $psi <= 180 && $psi >= -180 && $omega <= 180 && $omega >= -180) {
				$tmpphicol = $i;
				if (!defined $phicol) {
					$phicol = $tmpphicol;
					print "Assumming torsions start from column ".($i+1)."\n";
				}
				last;
			}
			if ((defined $phicol && $tmpphicol != $phicol) || !defined $tmpphicol) {
				close(NF);
				close(F);
				unlink($torsions_vall);
				die "ERROR! there was an issue when trying to determine the torsions columns from $vall\n";
			}
		}
		if (!defined $phicol) {
			close(NF);
			close(F);
			unlink($torsions_vall);
			die "ERROR! there was an issue when trying to determine the torsions columns from $vall\n";
		}
	}
	$vallline++;
	$vall{"$cols[0] $cols[1] $cols[2] $cols[3]"} = { line => $vallline, aa => $cols[1], ss => $cols[2], phi => $cols[$phicol], psi => $cols[$phicol+1],  omega => $cols[$phicol+2] };
}
close(F);

if ($opts{create_torsions_vall}) {
	my $torsions_vall = "$vall.torsions";
	print "Creating torsions vall: $torsions_vall\n";
	open(NF, ">$torsions_vall") or die "ERROR! cannot open $torsions_vall: $!\n";
	if ($vall =~ /(\.gz|\.Z)\s*$/) {
		open(F, "gzip -dc $vall |") or die "ERROR! cannot open $vall: $!\n";
	} else {
		open(F, $vall) or die "ERROR! cannot open $vall: $!\n";
	}
	while (<F>) {
		next if (/^#/);
		my @cols = split(/\s+/,$_);
		print NF $cols[1]." ".$cols[2]." ".$cols[$phicol]." ".$cols[$phicol+1]." ".$cols[$phicol+2]."\n";
	}
	close(F);
	close(NF);
}

foreach my $frag (@ARGV) {
	my ($pos, $frags, $phi, $psi, $omega, $code, $chain, $seqrespos, $aa, $ss, @frag);
	print "Creating indexed fragment: $frag.index\n";
	if ($frag =~ /(\.gz|\.Z)\s*$/) {
		open (F, "gzip -dc $frag |") or die "ERROR! cannot open $frag: $!\n";
	} else {
		open(F, $frag) or die "ERROR! cannot open $frag: $!\n";
	}
	open(NF, ">$frag.index") or die "ERROR! cannot open $frag.index: $!\n";
	print NF "# index 1 $vallline 0 $vallname\n";
	my $fragcnt = 0;
	while (<F>) {
		if (/^\s*position:\s+(\d+)\s+neighbors:\s+(\d+)\s*$/) {
			if ($fragcnt) {
				($fragcnt == $frags) or warn "WARNING! neighbors $frags does not equal the number of fragments $fragcnt for position $pos\n";
				if ($debug) { print "frags $fragcnt\n"; }
				print NF "\n";
			}
			$pos = $1;
			$frags = $2;
			$fragcnt = 0;
			next;
		} elsif (/^\s*$/ && $pos) {
			$code = undef;
			if (scalar@frag) {
				my $fragsize = scalar@frag;
				my $vallline = $vall{"$frag[0]->{code}$frag[0]->{chain} $frag[0]->{aa} $frag[0]->{ss} $frag[0]->{seqrespos}"}->{line};
				if ($debug) { print "fragsize $fragsize vallline: $vallline\n"; }
				my $output_torsions = 0;
				foreach my $f (@frag) {
					my $key = "$f->{code}$f->{chain} $f->{aa} $f->{ss} $f->{seqrespos}";
					(defined $vall{$key} && $vall{$key}->{phi} == $f->{phi} && $vall{$key}->{psi} == $f->{psi} && $vall{$key}->{omega} == $f->{omega} ) or do {
						warn "WARNING! there was an issue reading the fragment file $frag: no vall entry for $key, format error? Printing torsions.\n";
						$output_torsions = 1;
					};
					if ($debug) {
						print "#$key# fragsize: $fragsize $vall{$key}->{phi} == $f->{phi} && $vall{$key}->{psi} == $f->{psi} && $vall{$key}->{omega} == $f->{omega} vall line: $vall{$key}->{line}\n";
					}
				}
				($vallline) or do {
					warn "WARNING! there was an issue reading the fragment file $frag: no vall line entry for $frag[0]->{code}$frag[0]->{chain} $frag[0]->{aa} $frag[0]->{ss} $frag[0]->{seqrespos}, format error? Printing torsions.\n";
					$output_torsions = 1;
				};
				$fragcnt++;
				if ($output_torsions) {
					print NF "0 $fragsize\n";
					foreach my $f (@frag) {
						printf NF "$f->{aa} $f->{ss} $f->{phi} %9.3f %9.3f\n", $f->{psi}, $f->{omega};
					}
				} else {
					print NF "$vallline $fragsize\n";
				}
			}
			@frag = ();
		} elsif (!$code || scalar@frag) {
			if (/^\s*(\S{4})\s(\S+)\s+(\d+)\s(\w)\s(\w)\s+(\-?\d{1,3}\.\d{3})\s+(\-?\d{1,3}\.\d{3})\s+(\-?\d{1,3}\.\d{3})\s+/) {
				$code =  $1; $chain = $2; $seqrespos = $3; $aa = $4; $ss = $5; $phi = $6; $psi = $7; $omega = $8;
				push(@frag, { code => $code, chain => $chain, seqrespos => $seqrespos, aa => $aa, ss => $ss, phi => $phi, psi => $psi, omega => $omega });
			} else {
				close(F);
				close(NF);
				unlink("$frag.index");
				die "ERROR! there was an issue reading the fragment file $frag: frag format error?\n";
			}
		}
	}
	if ($fragcnt) {
		($fragcnt == $frags) or warn "WARNING! neighbors $frags does not equal the number of fragments $fragcnt for position $pos\n";
		if ($debug) { print "frags $fragcnt\n"; }
		print NF "\n";
	}
	close(F);
	close(NF);
}







sub getCommandLineOptions {
	use Getopt::Long;
	my $usage = qq{usage: $0
\t -vall		<vall file>
\t[-create_torsions_vall]

 <frag1> <frag2> ....

};

	# Get args
	my %opts = ();
	&GetOptions(
		\%opts,			 "vall=s",
		"create_torsions_vall"
	);

	# Check for legal invocation
	if (   !defined $opts{vall} || !scalar@ARGV )
	{
		print STDERR "$usage\n";
		exit(1);
	}

	# defaults
	#
	$opts{create_torsions_vall}		 ||= 0;

	return %opts;
}

__END__




rosetta++ fragment format examples:

 position:            1 neighbors:          200

 3erk _   252 S L  -69.350  175.247  178.303   -0.590    4.442    47.021 3     0.000 P  1 F  1
 3erk _   253 Q H  -62.828  -48.955  177.509   -0.590    4.442    47.021 3     0.000 P  1 F  1
 3erk _   254 E H  -50.940  -62.060 -179.282   -0.590    4.442    47.021 3     0.000 P  1 F  1

 1ehi A   195 N L -138.229  170.729 -176.893   -0.582    4.457    43.873 3     0.000 P  1 F  2
 1ehi A   196 A H  -68.569  -38.480 -179.597   -0.582    4.457    43.873 3     0.000 P  1 F  2
 1ehi A   197 E H  -65.993  -36.467  177.822   -0.582    4.457    43.873 3     0.000 P  1 F  2

 1b0n B    12 D L  -70.461  126.175 -173.552   -0.570    4.395     8.996 3     0.000 P  1 F  3
 1b0n B    13 Q H  -66.735  -27.413  175.997   -0.570    4.395     8.996 3     0.000 P  1 F  3
 1b0n B    14 E H  -71.444  -43.708  177.677   -0.570    4.395     8.996 3     0.000 P  1 F  3


position:            1 neighbors:          200

 3kda A   140 F L  -72.129  140.597 -175.348   14.250    8.200   30.120
 3kda A   141 P L  -75.033  151.596  175.608   14.690   11.960   29.840
 3kda A   142 A L  -84.425  -16.750 -177.429   13.010   14.060   27.160






