#!/usr/bin/perl
 
# extract the chains needed for docking: input the chains
# needed, with a dash between partners.

# JJG 1/3/2


# init - get the chains, make them upper case
$chains[0] = shift @ARGV;
$chains[0] = ($chains[0] eq '_') ? ' ' : uc $chains[0];
# chains[0] is all chains, chains[1] and [2] are the two partners

if (!$chains[0]) {
    print STDERR "usage: $0 <chains> [<pdb>]\n";
    print STDERR "  returns ATOM entries with the selected chains plus some header info\n";
    print STDERR "  chains should be given like AB-EF with a dash to indicate docking partners\n";
    print STDERR "  TERmini are placed after each partner only\n";
    exit -1;
}

# extract chains for partners 1 and 2
($tmp = $chains[0]) =~ /-/;
$chains[1] = $`;
$chains[2] = $';
$chains[1] = $chains[0] if !$chains[1] && !$chains[2];

# make regular expressions to search for all chains in the group  
for ($i=0; $i<3; ++$i) {
    ($chains_regexp[$i] = $chains[$i]) =~ s/([A-Z0-9\-])([A-Z0-9\-])/$1\|$2/g ;
    $chains_regexp[$i] =~ s/([A-Z0-9\-])([A-Z0-9\-])/$1\|$2/g ; # now have A|B|-|C|D
}

@pdb_buf = <>;

# header - keep lines beginning with HEADER COMPND TITLE and SOURCE
if ($chains[1]) {
    for ($i=0; $i <= 20; ++$i) {
	if ($pdb_buf[$i] =~ /^HEADER|^COMPND|^SOURCE|^TITLE/) {
	    print $pdb_buf[$i];
	}
    } # keep my comments too, after the header
    for ($i=0; $i <= 20; ++$i) {
	if ($pdb_buf[$i] =~ /^REMJJG/) {
	    print $pdb_buf[$i];
	}
    }
}

# detect chains in the file
%seen = (); $ch='';
for $line (@pdb_buf) {
    $ch=substr($line,21,1) if ($line =~ /^ATOM/);
    push @allchains,$ch unless $seen{$ch}++;
}

print STDERR "Extracting chains $chains[0] from @allchains\n";

print "REMJJG ---Extracted chains $chains[0] from @allchains\n";
print "REMJJG ---",scalar(localtime()),"\n"; 

# be sure desired chains are in the pdb
@desired_chains = split(//,$chains[0]);
$avail_chains = join '',@allchains,'-';
$avail_chains =~ s/([A-Z0-9\-])([A-Z0-9\-])/$1\|$2/g ;
$avail_chains =~ s/([A-Z0-9\-])([A-Z0-9\-])/$1\|$2/g ;
foreach $c (@desired_chains) {
    die("chain $c not found in @allchains") if ($c !~ $avail_chains);
}

# extract coordinates -- only ATOM and HETATM lines and last TER line
$body[1]=\@part1;
$body[2]=\@part2;
for ($p=1; $p<=2; ++$p) { # loop over each docking partner
    if ($chains[$p]) { # only collect the number of chains requested
	$lastter='TER';
	for ($i=0; $i <= $#pdb_buf; ++$i) { # loop through pdb
	    if ((substr ($pdb_buf[$i], 21, 1) =~ $chains_regexp[$p])) { # chain
		push @{ $body[$p] }, $pdb_buf[$i]
		    if ($pdb_buf[$i] =~ /^ATOM|^HETATM/);
		$lastter=$pdb_buf[$i]     
		    if ($pdb_buf[$i] =~ /^TER/);
	    }
	}
	chomp($lastter);
	push @{ $body[$p] },$lastter,"\n"; # only one terminus per partner
    }
}

print STDERR "First partner, $chains[1], $#part1 atoms\n";
print STDERR "Second partner, $chains[2], $#part2 atoms\n";

# smaller partner should be second in the file
if ($#part1 > $#part2 || ! $chains[1]) {
    print @part1,@part2;
}
else { # ($#part1 <= $#part2) { 
    print STDERR "Reversing order to $chains[2]-$chains[1]\n";
    print "REMJJG ---Reversing order to $chains[2]-$chains[1]\n";
    print @part2,@part1;
}

exit 0;

