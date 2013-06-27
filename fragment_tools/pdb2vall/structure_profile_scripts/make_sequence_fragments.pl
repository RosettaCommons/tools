#!/usr/bin/perl
use FindBin qw($Bin);

$| = 1;                                     # disable stdout buffering

my $VERSION = 1;

my $ROSETTA_PATH = "$Bin/../../../../";
my $picker_num_cpus = 8;

# read config
open(F, "$Bin/../pdb2vall.cfg") or die "ERROR! cannot open config file $Bin/../pdb2vall.cfg: $!\n";
while (<F>) {
  if (/^\s*rosetta_path\s*=\s*(\S+)\s*/) {
    $ROSETTA_PATH = $1;
  } elsif (/^\s*fragment_picker_num_cpus\s*=\s*(\d+)\s*/) {
    $picker_num_cpus = $1;
  }
}
close(F);
$ROSETTA_PATH =~ s/[\/\\]$//;
(-d $ROSETTA_PATH) or die "ERROR! rosetta installation path $ROSETTA_PATH does not exist\n";

###############################################################################
# config
###############################################################################

my $PICKER_NUM_CPUS   = $picker_num_cpus;
my $PICKER            = "$ROSETTA_PATH/main/source/bin/fragment_picker.boost_thread.linuxgccrelease";
my $ROSETTA_DATABASE  = "$ROSETTA_PATH/main/database";
my $VALL = "$Bin/residue_depth_vall/residue_depth_full_cullpdb_pc30_res2.0_R0.3_d110930_chains6280.4278.vall.gz";

###############################################################################
# init
###############################################################################

if (!-s $PICKER) {
	$PICKER =~ s/fragment_picker\.boost_thread/fragment_picker/;
}
(-s $PICKER) or die "ERROR! $PICKER does not exist.\n";
(-d $ROSETTA_DATABASE) or die "ERROR! $ROSETTA_DATABASE does not exist.\n";
(-s $VALL) or die "ERROR! $VALL does not exist.\n";

# argv
my %opts = &getCommandLineOptions ();
my $fasta = $opts{f};
my @fragsizes = ( 9 );

my $n_frags = 200;
my $n_candidates = 1000;

my $DEBUG = $opts{verbose};

my $native;
if (defined($opts{native})) {
  $native = $opts{native};
} else {
  die "ERROR! PDB file must be given using the -native option\n";
}

my $depth;
if (defined($opts{depth})) {
  $depth = $opts{depth};
} else {
  die "ERROR! residue depth file must be given using the -depth option\n";
}

if (defined($opts{frag_sizes})) {
  @fragsizes = split(/,/, $opts{frag_sizes});
  (!$DEBUG) || print "fragment sizes: ".join(" ", @fragsizes)."\n";
}
my $frags = join(" ", @fragsizes);

if (defined($opts{n_frags})) {
  $n_frags = $opts{n_frags};
  (!$DEBUG) || print "n_frags: $n_frags\n";
}

if (defined($opts{n_candidates})) {
  $n_candidates = $opts{n_candidates};
  (!$DEBUG) || print "n_candidates: $n_candidates\n";
}

(!$DEBUG) || print "FASTA: $fasta\n";

my $id;
if ($fasta =~ /(\w+)\.\w+\s*$/) {
  $id = $1;
}
if (!defined $id) {
  die "ERROR! id cannot be parsed from $fasta\n";
}
(!$DEBUG) || print "ID: $id\n";

###############################################################################
# main
###############################################################################

my $scorescfg =<<SCORESCFG;
# score name          priority    wght   min_allowed  extras
FragmentCrmsdResDepth      100     1.0       -
FragmentCrmsd               30     0.0       -
FragmentDME                 30     0.0       -
SCORESCFG

open(SCORESCFG_DEFS, ">picker_sequence_scores_$id.cfg");
print SCORESCFG_DEFS $scorescfg;
close(SCORESCFG_DEFS);

my $cmdtxt=<<CMDTXT;
-frags::write_sequence_only
-in:file:fasta		$fasta
-frags::describe_fragments     $id\_frags.fsc
-in::path::database     $ROSETTA_DATABASE
-in::file::vall         $VALL
-frags::n_candidates	$n_candidates
-frags::n_frags		$n_frags
-frags::frag_sizes	$frags
-out::file::frag_prefix $id
-frags::scoring::config picker_sequence_scores_$id.cfg
-frags::depth           $depth
-in:file:s              $native
CMDTXT

open(PATH_DEFS, ">picker_sequence_cmd_$id.txt");
print PATH_DEFS $cmdtxt;
close(PATH_DEFS);

my $time = time();
my $shell = "$PICKER \@picker_sequence_cmd_$id.txt -ignore_unrecognized_res true -j $PICKER_NUM_CPUS >& picker_sequence_cmd_$id.log";
print "shell: $shell\n";
system($shell);
my $diff = time() - $time;
print "RUNTIME: $diff seconds\n";

# done
exit 0;

###############################################################################
# util
###############################################################################

# getCommandLineOptions()
#
#  rets: \%opts  pointer to hash of kv pairs of command line options
#
sub getCommandLineOptions {
    use Getopt::Long;

    my $usage = qq{usage: $0
		   \t-verbose  specify for chatty output
                   \t-frag_sizes <size1,size2,...n>
                   \t-n_frags <number of fragments>
                   \t-n_candidates <number of candidates>
                   \t-depth <residue depth file>
                   \t-native <native pdb>
		   \t<fasta file>
		   \n\nVersion: $VERSION\n};

    # Get args
    my %opts = ();
    &GetOptions (\%opts, "verbose!","frag_sizes=s", "n_frags=i", "n_candidates=i", "depth=s", "native=s");


    if (scalar(@ARGV) != 1) {
      die "$usage\n";
    }

    $opts{f} = $ARGV[0];

    &checkExist("f",$opts{f});

    die("Fragment sizes are invalid\n") if (defined $opts{frag_sizes} && $opts{frag_sizes} !~ /^[\d\,]+$/);

    return %opts;
}

# checkExist()
#
sub checkExist {
    my ($type, $path) = @_;
    if ($type eq 'd') {
	if (! -d $path) {
            print STDERR "$0: dirnotfound: $path\n";
            exit -3;
	}
    }
    elsif ($type eq 'f') {
	if (! -f $path) {
            print STDERR "$0: filenotfound: $path\n";
            exit -3;
	}
	elsif (! -s $path) {
            print STDERR "$0: emptyfile: $path\n";
            exit -3;
	}
    }
}

