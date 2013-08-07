#!/usr/bin/perl

use Bio::AlignIO;
use FindBin qw($Bin);

# possible blastpgp locations
my @blastpgp = (
  "$Bin/../../../../../../src/blast/bin/blastpgp", # Robetta location
  "$Bin/../../blast/bin/blastpgp", # fragment_tools location
  "$Bin/../pdb_scripts/blastpgp"  # pdb2vall location
);

my @formatdb = (
	"$Bin/../../../../../../src/blast/bin/formatdb", # Robetta location
	"$Bin/../../blast/bin/formatdb", # fragment_tools location
	"$Bin/../pdb_scripts/formatdb"  # pdb2vall location
);

# find blastpgp
my $BLASTPGP;
foreach my $b (@blastpgp) {
  if (-s $b) {
    $BLASTPGP = $b;
    last;
  }
}

my $FORMATDB;
foreach my $f (@formatdb) {
	if (-s $f) {
		$FORMATDB = $f;
		last;
	}
}

# read config
open(F, "$Bin/../pdb2vall.cfg") or die "ERROR! cannot open config file $Bin/../pdb2vall.cfg: $!\n";
while (<F>) {
  if (/^\s*blastpgp\s*=\s*(\S+)\s*/) {
    $BLASTPGP = $1;
  }
	if (/^\s*formatdb\s*=\s*(\S+)\s*/) {
		$FORMATDB = $1;
	}
}
close(F);

(-s $BLASTPGP) or die "ERROR! blastpgp $BLASTPGP does not exist\n";

if (scalar@ARGV != 2) {
  print "USAGE: $0 <query fasta> <fasta alignment>\n";
  exit(1);
}
my ($queryfasta,$fasta) = @ARGV;
chomp $queryfasta;
chomp $fasta;


# Just a placeholder since we just need the binary checkpoint file that blastpgp initially generates
# given the input alignment.
my $placeholder = "$Bin/placeholder_seqs/placeholder_seqs";
if (!-s "$placeholder.phr") {
	system("$FORMATDB -o T -i $placeholder");
}

# get the sequence from the fasta file
my $sequence;
open(SEQFILE, "$queryfasta") or die "ERROR! cannot open $queryfasta: $!\n";
my $has_comment = 0;
my $has_eof     = 0;
my $eof;
my $header;
while (<SEQFILE>) {
  $eof = $_;
  if (/^>/) { $header = $_; $has_comment = 1; next; }
  s/\s//g;
  chomp;
  $sequence .= $_;
}
close(SEQFILE);
$has_eof = 1 if ($eof =~ /\n$/);
($has_comment && $has_eof) or die "fasta file must have a comment and end with a new line!\n";
$sequence =~ s/\s+//gs;

my $fastanew = $fasta.".new";
open(FA, $fasta) or die "ERROR! cannot open $fasta: $!\n";
open(FN, ">$fastanew") or die "ERROR! cannot open $fastanew: $!\n";
print FN $header;
print FN $sequence."\n";
foreach my $line (<FA>) {
  print FN $line;
}
close(FN);
close(FA);

my $blast = "$fastanew.blast";
my $tmpblast = "$fastanew.blast.tmp";

if (-s "$blast.checkpoint") {
  print "$blast.checkpoint exists\n";
  print "done!\n";
  exit(0);
}


my $in  = Bio::AlignIO->new(-file   => "$fastanew", -format => 'fasta');
my $out = Bio::AlignIO->new(-file   => ">$tmpblast" ,  -format => 'psi');
while ( my $aln = $in->next_aln() ) { $out->write_aln($aln); }
open(F, "$blast.tmp");
open(NF, ">$blast");
while (<F>) {
  $_ =~ s/\./-/gs; # need to switch . to '-' for blastpgp
  print NF $_;
}
close(F);
close(NF);
unlink("$blast.tmp");

# generate blast binary checkpoint for PSSM from alignment
my $check = "$blast.chk";
my $blastcmd = "( $BLASTPGP -i $queryfasta -B $blast -t 1 -j 1 -h 0.001 -e 0.001 -b 0 -k 0 -C $check -d $placeholder  ) >& /dev/null";
print "$blastcmd\n";
(system($blastcmd) == 0 && -s $check) or die "ERROR! blast failed\n";


# parse the checkpoint matrix.
my @checkpoint_matrix = &parse_checkpoint_file("$check");
@checkpoint_matrix = &finish_checkpoint_matrix($sequence, @checkpoint_matrix);
&write_checkpoint_file("$blast.checkpoint", $sequence, @checkpoint_matrix);

print "done!\n";

## parse_checkpoint_file -- parses a PSI-BLAST binary checkpoint file.
#
# args:  filename of checkpoint file.
# rets:  N x 20 array containing checkpoint weight values, where N
#        is the size of the protein that BLAST thought it saw...

sub parse_checkpoint_file {
  my $filename = shift;
  my $buf;
  my $seqlen;
  my $seqstr;
  my $i;
  my $j;
  my @aa_order = split(//,'ACDEFGHIKLMNPQRSTVWY');
  my @altschul_mapping = (0,4,3,6,13,7,8,9,11,10,12,2,14,5,1,15,16,19,17,18);
  my @w;
  my @output;

  open(INPUT, $filename) || die ("Couldn't open $filename for reading.\n");

  read(INPUT, $buf, 4) || die ("Couldn't read $filename!\n");
  $seqlen = unpack("i", $buf);

  (!$DEBUG) || print "Sequence length: $seqlen\n";

  read(INPUT, $buf, $seqlen) || die ("Premature end: $filename.\n");
  $seqstr = unpack("a$seqlen", $buf);

  for ($i = 0; $i < $seqlen; ++$i) {
    read(INPUT, $buf, 160) || die("Premature end: $filename, line: $i\n");
    @w = unpack("d20", $buf);

    for ($j = 0; $j < 20; ++$j) {
      $output[$i][$j] = $w[$altschul_mapping[$j]];
    }
  }

  return @output;
}

## finish_checkpoint_matrix -- "finishes" the parsed PSI-BLAST checkpoint matrix,
##                             by adding pseudo-counts to any empty columns.
#
# args:  1) sequence string  2) array returned by parse_checkpoint_file
# rets:  "finished" array.  suitable for printing, framing, etc.

sub finish_checkpoint_matrix {
  my ($s, @matrix) = @_;
  my @sequence = split(//,$s);
  my $i;
  my $j;
  my $line;
  my $sum;
  my @words;
  my @b62;
  my @blos_aa = (0,14,11,2,1,13,3,5,6,7,9,8,10,4,12,15,16,18,19,17);
  my %aaNum = ( A => 0,
                C => 1,
                D => 2,
                E => 3,
                F => 4,
                G => 5,
                H => 6,
                I => 7,
                K => 8,
                L => 9,
                M => 10,
                N => 11,
                P => 12,
                Q => 13,
                R => 14,
                S => 15,
                T => 16,
                V => 17,
                W => 18,
                Y => 19,
                X => 0   ### cheep fix for now ##
              );

  (length($s) == scalar(@matrix)) || die ("Length mismatch between sequence and checkpoint file!\n");

  my $blosum62 = <<BLOSUM;
#  BLOSUM Clustered Target Frequencies=qij
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
   A      R      N      D      C      Q      E      G      H      I      L      K      M      F      P      S      T      W      Y      V
0.0215
0.0023 0.0178
0.0019 0.0020 0.0141
0.0022 0.0016 0.0037 0.0213
0.0016 0.0004 0.0004 0.0004 0.0119
0.0019 0.0025 0.0015 0.0016 0.0003 0.0073
0.0030 0.0027 0.0022 0.0049 0.0004 0.0035 0.0161
0.0058 0.0017 0.0029 0.0025 0.0008 0.0014 0.0019 0.0378
0.0011 0.0012 0.0014 0.0010 0.0002 0.0010 0.0014 0.0010 0.0093
0.0032 0.0012 0.0010 0.0012 0.0011 0.0009 0.0012 0.0014 0.0006 0.0184
0.0044 0.0024 0.0014 0.0015 0.0016 0.0016 0.0020 0.0021 0.0010 0.0114 0.0371
0.0033 0.0062 0.0024 0.0024 0.0005 0.0031 0.0041 0.0025 0.0012 0.0016 0.0025 0.0161
0.0013 0.0008 0.0005 0.0005 0.0004 0.0007 0.0007 0.0007 0.0004 0.0025 0.0049 0.0009 0.0040
0.0016 0.0009 0.0008 0.0008 0.0005 0.0005 0.0009 0.0012 0.0008 0.0030 0.0054 0.0009 0.0012 0.0183
0.0022 0.0010 0.0009 0.0012 0.0004 0.0008 0.0014 0.0014 0.0005 0.0010 0.0014 0.0016 0.0004 0.0005 0.0191
0.0063 0.0023 0.0031 0.0028 0.0010 0.0019 0.0030 0.0038 0.0011 0.0017 0.0024 0.0031 0.0009 0.0012 0.0017 0.0126
0.0037 0.0018 0.0022 0.0019 0.0009 0.0014 0.0020 0.0022 0.0007 0.0027 0.0033 0.0023 0.0010 0.0012 0.0014 0.0047 0.0125
0.0004 0.0003 0.0002 0.0002 0.0001 0.0002 0.0003 0.0004 0.0002 0.0004 0.0007 0.0003 0.0002 0.0008 0.0001 0.0003 0.0003 0.0065
0.0013 0.0009 0.0007 0.0006 0.0003 0.0007 0.0009 0.0008 0.0015 0.0014 0.0022 0.0010 0.0006 0.0042 0.0005 0.0010 0.0009 0.0009 0.0102
0.0051 0.0016 0.0012 0.0013 0.0014 0.0012 0.0017 0.0018 0.0006 0.0120 0.0095 0.0019 0.0023 0.0026 0.0012 0.0024 0.0036 0.0004 0.0015 0.0196
BLOSUM

  $i = 0;
  my @blosum62 = split( /\n/, $blosum62 );

  # read/build the blosum matrix
  foreach my $blosumline (@blosum62) {
    next if ($blosumline !~ /^\d/);
    chomp;
    @words = split(/\s/, $blosumline);

    for ($j = 0; $j <= $#words; ++$j) {
      $b62[$blos_aa[$i]][$blos_aa[$j]] = $words[$j];
      $b62[$blos_aa[$j]][$blos_aa[$i]] = $words[$j];
    }

    ++$i;
  }

  # normalize the blosum matrix so that each row sums to 1.0
  for ($i = 0; $i < 20; ++$i) {
    $sum = 0.0;

    for ($j = 0; $j < 20; ++$j) {
      $sum += $b62[$i][$j];
    }

    for ($j = 0; $j < 20; ++$j) {
      $b62[$i][$j] = ($b62[$i][$j] / $sum);
    }
  }

  # substitute appropriate blosum row for 0 rows
  for ($i = 0; $i <= $#matrix; ++$i) {
    $sum = 0;

    for ($j = 0; $j < 20; ++$j) {
      $sum += $matrix[$i][$j];
    }

    if ($sum == 0) {
      for ($j = 0; $j < 20; ++$j) {
        $matrix[$i][$j] = $b62[$aaNum{$sequence[$i]}][$j];
      }
    }
  }

  return @matrix;
}

sub write_checkpoint_file {
  my ($filename, $s, @matrix) = @_;
  my $row;
  my $col;
  my @seq = split(//,$s);

  open(OUTPUT, ">$filename");

  die ("Length mismatch between sequence and checkpoint matrix!\n") unless (length($s) == scalar(@matrix));

  print OUTPUT scalar(@matrix),"\n";

  for ($row = 0; $row <= $#matrix; ++$row) {
    print OUTPUT "$seq[$row] ";
    for ($col = 0; $col < 20; ++$col) {
      printf OUTPUT "%6.4f ", $matrix[$row][$col];
    }
    print OUTPUT "\n";
  }

  print OUTPUT "END";

  close(OUTPUT);
}



