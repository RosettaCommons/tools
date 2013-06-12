#!/usr/bin/perl
# determines for a run for each cluster the pdb file with the lowest score:
# input: clusteroutput of rms table (from modified chris script rms2avglink.csh) X_TAR_top200_Y.rms
#        score file  copied from /data/morozov/CAPRI/X_TAR/top200_Y_scores
# output: list of ranks for each pdbfile in cluster

die ("Usage: $0 <clusterfile> <scorefile>") if (@ARGV<1);
$clusterfile = $ARGV[0];
$scorefile = $ARGV[1];

# score to use:
# 1,3: JG 2,4,5:TK - TK=2 JG=41 for target 2,3 JG=9 for target 1
$score_column = 2;
$rms_column = 3;
if ($clusterfile =~ /_TAR_top200_1\.cluster/ || $clusterfile =~ /_TAR_top200_3\.cluster/){
    if ($clusterfile =~ /1_TAR/){
	$score_column = 9;
    }else{
	$score_column = 41;
    }
}else { if ($scorefile =~ /reweight/) {$score_column = 34; }}


printf STDERR "score_column=$score_column\n";
printf STDERR "rms_column=$rms_column\n";

# read clusterfile
open (SCORES,"$scorefile")||die ("scorefile $scorefile could not be opened");
while ($line = <SCORES>){
    chomp($line);
    @words = split /\s+/,$line;
    $name = $words[0];
# strip name of run (after .pdb) from name
# if same name: comment out
    @tmp = split(/pdb/,$name);
    $name = $tmp[0]."pdb";
#
    $scores{$name}  = $words[$score_column-1];
    $rmss{$name}  = $words[$rms_column-1];
#    print "$name $scores{$name}\n";
}
close (SCORES);

print "#cluster number -  cluster size - best decoy - best score - rms\n";

# read clusterfile
open (CLUSTERS,"$clusterfile")||die ("file $clusterfile could bot be opened");
$previous_index = 0;
$size = 0;
$minimum = 9999999;
$best = "none"; 

while ($line = <CLUSTERS>){
  chomp($line);
  next if ($line =~ /pdbs/); # skip header if present
  print "attention!! $line\n" if ($line !~ /^(\w+\.\w+\.pdb)\s+(\d+)/);
  $pdbname = $1;
  $index = $2;
  $score = $scores{$pdbname};
  $rms = $rmss{$pdbname};
  if (!defined($scores{$pdbname})){
      print "attention: no score available $pdbname: skipped $pdbname\n";
      next;
  }
  if ($index != $previous_index){
      printf "%3d %3d %20s %7.2f %6.2f\n",$previous_index,$size,$best,$minimum,$bestrms if ($previous_index>0);
      $previous_index = $index;
      $size = 0;
      $minimum = 9999999;
      $best = "none"; 
  }
  $size++;
  if ($score<$minimum){
      $minimum = $score;
      $best = $pdbname;
      $bestrms = $rms;
  }
}
close(CLUSTERS);
printf "%3d %3d %20s %7.2f %6.2f\n",$previous_index,$size,$best,$minimum,$bestrms;









