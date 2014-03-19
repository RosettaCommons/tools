#!/usr/bin/perl -D

use CGI;

	$parse = new CGI;
	
	$sequence = $parse->param('sequence');
	$num = $parse->param('num');

	if($sequence) {
		print_output();
	}
	else {
		print_form();
	}
	
	exit;

####################################################
## Begin Subs
###################################################

### form sub ###
sub print_form {

    print <<EOF;
Content-type: text/html

<html><head><title>Data Input</title></head>
<body bgcolor=lightgrey>
<b>Please enter multiple sequence alignment (MSA) of TM helix in Fasta format:</b><br>
<form method=POST action="$ENV{SCRIPT_NAME}">
<textarea cols=40 rows=20 name=sequence></textarea>
<br><br><br>
<b>Number of the first residue in MSA:</b><br>
<input type="text" name=num>
<input type=submit value=submit>
</form>

</body></html>
EOF
}

#################################################################
### main output sub ###
#################################################################

sub print_output {

print <<EOF;
Content-type: text/html

<html><head><title>Results</title></head>
<body bgcolor=lightgrey>
<b>Results:</b><br><pre>
EOF


## Changes made by mdabro3 - refer to the entropy_W.pl for original script

### not sure what to do here
#$alignment = $ARGV[0];

#variable "resnum" here is a sequence number of the first residue of the query sequence in a seq alignment
### puts user input into the place we are looking for it
$resnum = $num;

$n=0;
$sump=0;
$sumf=0;
$sumim=0;
$aanum=0;

@amino = ( "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K",
  "M", "F", "P", "S", "T", "W", "Y", "V" 
);

%propi = (
        A => 0.71,
        R => 1.47,
        N => 0.96,
        D => 1.20,
        C => 1.16,
        Q => 0.61,
        E => 0.90,
        G => 0.48,
        H => 0.82,
        I => 1.11,
        L => 1.18,
        K => 2.38,
        M => 1.38,
        F => 1.57,
        P => 0.99,
        S => 0.69,
        T => 0.72,
        W => 2.45,
        Y => 1.23,
        V => 0.98 
);

%propm = (
        A => 0.82,
        R => 0.18,
        N => 0.19,
        D => 0.29,
        C => 1.01,
        Q => 0.26,
        E => 0.19,
        G => 0.35,
        H => 0.12,
        I => 1.88,
        L => 1.71,
        K => 0.42,
        M => 1.02,
        F => 1.97,
        P => 0.65,
        S => 0.55,
        T => 0.66,
        W => 1.65,
        Y => 0.94,
        V => 1.77
);


#calculate entropy values for each position of the alignment;

#open(DATA, $alignment) || die "couldn't open \"$data_name\": $!\n";
#while(<DATA>){
#chomp;
#	push @LoL, [ split// ];
#};

#lines added by mdabro3
@tmp =  split(/\s/, $sequence);

#debug info:
#$rows=$#tmp;
#print "lines after first split: $rows\n";

foreach (@tmp){
chomp;
	if ( /\w/ ){
		push @LoL, [ split// ];
	}
}

$nrow=$#LoL+1;
$len=$#{$LoL[$i]}+1;
$bnum = int($len/5);
#debug info:
#print "nrow=$nrow length=$len\n";
for $i ( 0 .. $#LoL ){
	for $j ( 0 .. $#{$LoL[$i]} ) {
		$r=$LoL[$i][$j];
		$oc{$r." ".$j}++;
	};
};

for $j ( 0 .. $#{$LoL[$i]} ) {
   foreach $a (@amino) {
      $p{$a}=$oc{$a." ".$j}/$nrow;
         if ($p{$a} ne 0){
            $h{$j}+=$p{$a}*log $p{$a};
            if(($j <= $bnum)||($j > $len-$bnum)){
               $lip{$j}+=$p{$a}*$propi{$a};
            }else{
               $lip{$j}+=$p{$a}*$propm{$a};
            };
         };
   };
	$ph{$j}=2.718**((-1)*$h{$j});
};
foreach $j ( sort keys(%ph)){
	$r=$LoL[0][$j];
	$m=$resnum+$j;
	$sump +=$ph{$j};
        $sumLIP +=$lip{$j};
};

$bnum = int($len/5);

# calculate average propensity and entropy values for each of 7 surfaces;
#Start with the first residue of the query sequence and get propensity and entropy for every 7th residue:

for($i=0; $i<4; $i++){
   print "     SURFACE $i:\n";
      for ($j=$i; $j<$len; $j +=7){
         $r=$LoL[0][$j];
         $sumim{$i} +=$lip{$j};
         $prop=$lip{$j};
         $sume{$i} += $ph{$j};
         $aanum{$i} ++;
         $rn=$j+$resnum;
         $r3 = & residuename123($r);

#print output:
printf "%3s %s %6.3f %5.3f\n", $rn, $LoL[0][$j], $prop, $ph{$j};

#Sample every 3rd and 4th residues between the 7th residues:

                for ($k=$j+3; $k<=$j+4; $k++){
                   if($k < $len){
                      $r=$LoL[0][$k];
                      $r3 = & residuename123($r);
                      $sumim{$i} +=$lip{$k};
                      $prop=$lip{$k};
                      $sume{$i} += $ph{$k};
                      $aanum{$i} ++;
                      $rn=$k+$resnum;

#print output:
printf "%3s %s %6.3f %5.3f\n", $rn, $LoL[0][$k], $prop, $ph{$k};
                   };
	       };
      };
};

for($i=4; $i<=6; $i++){
print "     SURFACE $i:\n";
        for ($j=$i; $j<$len; $j +=7){
                        $r=$LoL[0][$j];
                        $sumim{$i} +=$lip{$j};
                        $prop=$lip{$j};
                        $sume{$i} += $ph{$j};
                        $aanum{$i} ++;
                        $rn=$j+$resnum;
                        $r3 = & residuename123($r);

printf "%3s %s %6.3f %5.3f\n", $rn, $LoL[0][$j], $prop, $ph{$j};
#                   };

        for ($k=$j+3; $k<=$j+4; $k++){
           if($k < $len){
              $r=$LoL[0][$k];
              $r3 = & residuename123($r);
              $rn=$k+$resnum;
              $sumim{$i} +=$lip{$k};
              $prop=$lip{$k};
              $sume{$i} += $ph{$k};
              $aanum{$i} ++;

printf "%3s %s %6.3f %5.3f\n", $rn, $LoL[0][$k], $prop, $ph{$k};
            };
        };
      };
          for ($k=$i-4; $k<=$i-3; $k++){
             if($k < $len){
                $r=$LoL[0][$k];
                $r3 = & residuename123($r);
                $sumim{$i} +=$lip{$k};
                $prop=$lip{$k};
                $sume{$i} += $ph{$k};
                $aanum{$i} ++;
                $rn=$k+$resnum;

printf "%3s %s %6.3f %5.3f\n", $rn, $LoL[0][$k], $prop, $ph{$k};
              };
	   };
};
#print output result;
#print "$alignment\n";

print "<br><br>";
print "SURFACE LIPOPHILICITY ENTROPY   LIPS\n";

foreach $i (keys(%sumim)){
        $avpim=$sumim{$i}/$aanum{$i};
        $avpim=$avpim*2;
	$ave=$sume{$i}/$aanum{$i};
	$peim=$avpim*$ave;
        printf "%4s %12.3f %10.3f %8.3f\n", $i, $avpim, $ave, $peim;
};

#printf "<br><br><br>end</pre>";
printf "</body></html>";

}
####################
### end main

###################
### more subs


sub residuename321 {
  my($sa) = @_;
  my($res);
 if($sa eq "ALA") {
    $res = "A";
  }elsif($sa eq "ARG") {
    $res = "R";
  }elsif($sa eq "ASN") {
    $res = "N";
  }elsif($sa eq "ASP") {
    $res = "D";
  }elsif($sa eq "CYS") {
    $res = "C";
  }elsif($sa eq "GLY") {
    $res = "G";
  }elsif($sa eq "GLU") {
    $res = "E";
  }elsif($sa eq "GLN") {
    $res = "Q";
  }elsif($sa eq "HIS") {
    $res = "H";
  }elsif($sa eq "ILE") {
    $res = "I";
  }elsif($sa eq "LEU") {
    $res = "L";
  }elsif($sa eq "LYS") {
    $res = "K";
  }elsif($sa eq "MET") {
    $res = "M";
  }elsif($sa eq "PHE") {
    $res = "F";
  }elsif($sa eq "PRO") {
    $res = "P";
  }elsif($sa eq "SER") {
    $res = "S";
  }elsif($sa eq "THR") {
    $res = "T";
  }elsif($sa eq "TRP") {
    $res = "W";
  }elsif($sa eq "TYR") {
    $res = "Y";
  }else {
    $res = "V";
  }
  return $res;
}

sub residuename123 {
  my($sa) = @_;
  my($res);
 if($sa eq "A") {
    $res = "ALA";
  }elsif($sa eq "R") {
    $res = "ARG";
  }elsif($sa eq "N") {
    $res = "ASN";
  }elsif($sa eq "D") {
    $res = "ASP";
  }elsif($sa eq "C") {
    $res = "CYS";
  }elsif($sa eq "G") {
    $res = "GLY";
  }elsif($sa eq "E") {
    $res = "GLU";
  }elsif($sa eq "Q") {
    $res = "GLN";
  }elsif($sa eq "H") {
    $res = "HIS";
  }elsif($sa eq "I") {
    $res = "ILE";
  }elsif($sa eq "L") {
    $res = "LEU";
  }elsif($sa eq "K") {
    $res = "LYS";
  }elsif($sa eq "M") {
    $res = "MET";
  }elsif($sa eq "F") {
    $res = "PHE";
  }elsif($sa eq "P") {
    $res = "PRO";
  }elsif($sa eq "S") {
    $res = "SER";
  }elsif($sa eq "T") {
    $res = "THR";
  }elsif($sa eq "W") {
    $res = "TRP";
  }elsif($sa eq "Y") {
    $res = "TYR";
  }else {
    $res = "VAL";
  }
  return $res;
}


