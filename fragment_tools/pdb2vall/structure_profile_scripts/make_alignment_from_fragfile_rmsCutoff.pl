#!/usr/bin/perl

# for max fragment rank considered
my $TOPN = 25;

# for fragment quality
my $MAXCRMSD = 0.5;

# for chainbreaks ca-ca dist squared max 16 (4 angstroms)
my $MAXCADISTSQ = 16;

if (scalar@ARGV != 1) {
  print "USAGE: $0 <fasta (codechain.fasta)>\n"; # <sequence only fragfile> <frag score file>\n";
  exit(1);
}

my ($fasta) = @ARGV;
chomp $fasta;
my $frag = $fasta;
my $scorefile = $fasta;

$frag =~ s/([^\/]+)\.fasta/$1\.\*mers/;
$scorefile =~ s/([^\/]+)\.fasta/$1\_frags.fsc.\*mers/;

my @frags = `ls -1 $frag`;
my @scorefiles = `ls -1 $scorefile`;
if (scalar@frags !=1) { die "ERROR! cannot find $frag\n"; }
if (scalar@scorefiles !=1) { die "ERROR! cannot find $scorefile\n"; }

$frag = $frags[0];
$scorefile = $scorefiles[0];
chomp $frag;
chomp $scorefile;
(-s $frag) or die "ERROR! fragment file $frag must exist\n";
(-s $scorefile) or die "ERROR! fragment score file $scorefile must exist\n";

my $pdb = $fasta;
$pdb =~ s/\.fasta$/\.pdb/;
(-s $pdb) or die "ERROR! native PDB $pdb must exist\n";

my $finalout = $frag.".ali.fasta";
#if (-s $finalout) {
#  print "$finalout exists!\n";
#  exit(1);
#}

my $seq = "";
open(F, $fasta);
while (<F>) {
 next if (/^>/);
 $seq .= $_;
}
close(F);
$seq =~ s/\s+//gs;

#print "$seq\n";

my @seq = split(//, $seq);

#find chain breaks
open(F, $pdb);
my $res_i = 0;
my $info = [];
my %ca;
while (my $l = <F>) {
    if ($l =~ /^ATOM/) {
      my $name = substr($l, 12, 4);
      if ($name =~ /CA/) {
        my $aa = &three_to_one( substr($l,17,3) );
        my $chain = substr($l,21,1);
        my $pos = &trim(substr($l,22,4));
        if ( $chain !~ /^\s*$/ ) {
          $info->[$res_i]->{chain} = $chain;
        } else {
          $info->[$res_i]->{chain} = '_';
        }
	if (!$ca{$pos}) {
          $info->[$res_i]->{CA_x}  =  &trim(substr($l,30,8));
          $info->[$res_i]->{CA_y}  =  &trim(substr($l,38,8));
          $info->[$res_i]->{CA_z}  =  &trim(substr($l,46,8));
          $info->[$res_i]->{pos} = $pos;
          $info->[$res_i]->{aa} = $aa;
          ++$res_i;
	  $ca{$pos} = 1;
	}
      }
    }
}
close(F);

#query_pos  vall_pos  pdbid c ss  FragmentCrmsdResDepth FragmentCrmsd FragmentDME  TOTAL  FRAG_ID
#          1        210  1r89 A L                  0.95          0.66        0.42    0.952    963123

open(F, $scorefile) or die "ERROR! cannot open $scorefile: $!\n";;
while (my $l = <F>) {
  if ($l =~  /^\s*([\d\-]+)\s+/) {
    $l =~ s/^\s+//;
    my @cols = split(/\s+/,$l);
    my $qpos = $cols[0];
    my $crmsd = $cols[6];
    push(@{$qpos.'_crmsd'}, $crmsd);
  }
}
close(F);


my %chainbreak;
for (my $i=0;$i<$#seq;$i++) {
#  print "$i $info->[$i]->{aa} $seq[$i]\n";
  ($info->[$i]->{aa} eq $seq[$i]) or die "ERROR! aa does not match at pos ".($i+1)."  $info->[$i]->{aa} ne $seq[$i]\n";
  my $dist = ($info->[$i]->{CA_x}-$info->[$i+1]->{CA_x})**2 + ($info->[$i]->{CA_y}-$info->[$i+1]->{CA_y})**2 + ($info->[$i]->{CA_z}-$info->[$i+1]->{CA_z})**2;
  if ($dist > $MAXCADISTSQ) {
    #chain break
    $chainbreak{$i+1} = 1;
    warn "Warning! chain break at ".($i+1)." $dist\n";
   }
}




open(F, $frag);
my $fraglen = 0;
my $pos;
my $alicnt;
while (<F>) {
  if (/^position:\s+(\d+)\s+/) {
    $pos = $1;
    $alicnt = 0;
  } elsif (/^\s+(\w\w\w\w)\s+(\S)\s+(\d+)\s+(\w+)\s*$/) {
    $alicnt++;
    next if ($alicnt > $TOPN);
    my $pdb = $1; my $pdbchain = $2; my $startpos = $3;
    my $ali = $4;
    $fraglen = length($ali) if (!$fraglen);

    my $skip = 0;

    # skip fragments spanning a chainbreak
    foreach my $fpos ($pos .. $pos + $fraglen - 1) {
      if ($fpos < $pos + $fraglen - 1 && exists $chainbreak{$fpos} && $chainbreak{$fpos}) {
        $skip = 1; last;
      }
    }
    if (!$skip) {
      my @qposcrmsd = @{$pos.'_crmsd'};
      if ($qposcrmsd[$alicnt-1] > $MAXCRMSD) {
         warn "Warning! skipping fragment at pos $pos rank $alicnt: crmsd $qposcrmsd[$alicnt-1] > $MAXCRMSD\n";
         $skip = 1;
      }
    }

    if (!$skip) {
      my $finalali = '-'x($pos-1).$ali.('-'x(length($seq)-length($ali)-$pos+1));
      my @tmp = split(//, $finalali);
      for (my $i=0;$i<=$#tmp;$i++) {
        push( @{$i}, $tmp[$i] ) if ($tmp[$i] ne '-');
      }
#    print "$finalali\n";
    }
  }
}
close(F);
my $cnt = 1;
open(FO, ">$finalout");
foreach my $index (0 .. $fraglen*$TOPN+50) {
  my $final = "";
  for (my $i=0;$i<=$#seq;$i++) {
    my @col = @{$i};
    if ($#col >= $index) {
      $final .= $col[$index];
    } else { $final .= "-"; }
  }
  if ($final !~ /^\-+$/) {
    printf FO ">seq%4.4d\n", $cnt;
    print FO "$final\n";
    $cnt++;
  }
}
close(FO);

sub trim {
  my $str = shift;
  $str =~ s/^\s+|\s+$//g;
  return $str;
}

sub three_to_one {
  my $three_letter = shift;
  my %three_to_one = (
              'ALA' => 'A', 'ARG' => 'R', 'ASN' => 'N', 'ASP' => 'D',
              'CYS' => 'C', 'GLU' => 'E', 'GLN' => 'Q', 'GLY' => 'G',
              'HIS' => 'H', 'ILE' => 'I', 'LEU' => 'L', 'LYS' => 'K',
              'MET' => 'M', 'PHE' => 'F', 'PRO' => 'P', 'SER' => 'S',
              'THR' => 'T', 'TRP' => 'W', 'TYR' => 'Y', 'VAL' => 'V',
              'UNK' => 'X' );
  if (exists $three_to_one{$three_letter}) {
    return $three_to_one{$three_letter};
  } else {
    warn "WARNING! no three to one conversion for '$three_letter': using 'X'\n";
    return 'X';
  }
}


__END__

position:            1 neighbors:           50

 3aaf A   956 SWDFGPQAF
 3jq0 A   410 GEDPTEEIN
 2vbk A   172 VTDNYQAIQ
 1s1d A   189 HENWVSNYN
 3ckc A   453 GGDATGDIN
 3fdh A   371 NRDAEDALK
 2w61 A   112 TKSHDICME
 3ckc A   239 GQTDYAKAE
 3dhu A   332 HGDVTPLIQ
 3fgr B   415 VADKTAELY
 3p5p A   670 GRDMLAHIR

