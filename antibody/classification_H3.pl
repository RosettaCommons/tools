#!/bin/perl

use strict;
use warnings;

if(@ARGV != 2){
	print "USAGE: perl *.pl [LIST THAT CONTAINS ANTIBODY NAME, LENGTH OF H3 and DERECTORY WHERE YOUR DECOYS EXISTS] [NUMBER OF DECOYS]\n";
	exit;
}

my $list = $ARGV[0];
my $n_conf = $ARGV[1];
my $path;
my @data;
my $target;
my $len_h3;

my $num;
my $atom;
my $cid;
my $resnum;
my $x;
my $y;
my $z;
my $base;
my $dist = 0;
my $pre_dist;
my $x_o;
my $x_ne1;
my $y_o;
my $y_ne1;
my $z_o;
my $z_ne1;
my $cnt_k = 0;
my $cnt_e = 0;
my $cnt_o = 0;
my $cnt = 0;
my $base_start = 0;


open(RESULT, ">> ./classification_result");
print RESULT "#Antibody\tH3_length\tN of kink\tN of extend\tN of others\n";
open(LIST, $list) or die "check\t$list\n";
while(my $line = <LIST>){
	if($line =~ /^#/){
		next;
	}

	$line =~ s/\n//;
	@data = split(/\s+/, $line);
	$target  = $data[0];
	$len_h3  = $data[1];
	$path    = $data[2];

	if(-e "$path"){

	}else{
		print "No directory...$path\n";
		next;
	}

	if($len_h3 == 34){
		$base_start = "100Z";
	}elsif($len_h3 == 33){
		$base_start = "100Y";
	}elsif($len_h3 == 32){
		$base_start = "100X";
	}elsif($len_h3 == 31){
		$base_start = "100W";
	}elsif($len_h3 == 30){
		$base_start = "100V";
	}elsif($len_h3 == 29){
		$base_start = "100U";
	}elsif($len_h3 == 28){
		$base_start = "100T";
	}elsif($len_h3 == 27){
		$base_start = "100S";
	}elsif($len_h3 == 26){
		$base_start = "100R";
	}elsif($len_h3 == 25){
		$base_start = "100Q";
	}elsif($len_h3 == 24){
		$base_start = "100P";
	}elsif($len_h3 == 23){
		$base_start = "100O";
	}elsif($len_h3 == 22){
		$base_start = "100N";
	}elsif($len_h3 == 21){
		$base_start = "100M";
	}elsif($len_h3 == 20){
		$base_start = "100L";
	}elsif($len_h3 == 19){
		$base_start = "100K";
	}elsif($len_h3 == 18){
		$base_start = "100J";
	}elsif($len_h3 == 17){
		$base_start = "100I";
	}elsif($len_h3 == 16){
		$base_start = "100H";
	}elsif($len_h3 == 15){
		$base_start = "100G";
	}elsif($len_h3 == 14){
		$base_start = "100F";
	}elsif($len_h3 == 13){
		$base_start = "100E";
	}elsif($len_h3 == 12){
		$base_start = "100D";
	}elsif($len_h3 == 11){
		$base_start = "100C";
	}elsif($len_h3 == 10){
		$base_start = "100B";
	}elsif($len_h3 == 9){
		$base_start = "100A";
	}elsif($len_h3 == 8){
		$base_start = "100";
	}elsif($len_h3 == 7){
		$base_start = "99";
	}elsif($len_h3 == 6){
		$base_start = "98";
	}elsif($len_h3 == 5){
		$base_start = "97";
	}elsif($len_h3 == 4){
		$base_start = "96";
	}elsif($len_h3 == 3){
		$base_start = "95";
	}

	if(-e "out_classification/dist"){

	}else{
		`mkdir -p out_classification/dist`;
		`mkdir -p out_classification/angle`;
		`mkdir -p out_classification/base`;
	}


	if(-e "out_classification/base/$target.base"){
		`rm out_classification/base/$target.base`;
	}

	if(-e "out_classification/angle/$target\_0001.base"){
		`rm out_classification/angle/*`;
	}

	if(-e "out_classification/dist/$target\_0001.dist"){
		`rm out_classification/dist/*`;
	}

	print "$target\t$len_h3\t$path\n";

	open(ALL, ">out_classification/classify_$target.txt");

	for(my $i = 1; $i <= $n_conf; $i++){
		if($i < 10){
			$num = "000$i";
		}elsif($i < 100){
			$num = "00$i";
		}elsif($i < 1000){
			$num = "0$i";
		}elsif($i < 10000){
			$num = "$i";
		}

		if(-e "$path/$target\_$num.pdb" ){
			open(DIST, ">./out_classification/dist/$target\_$num.dist") or die "check2\n";
			open(DATA, ">./out_classification/angle/$target\_$num.base") or die "check3\n";
			open(PDB, "$path/$target\_$num.pdb") or die "check5 $target\t$num.pdb\n";
			while(my $line2 = <PDB>){
				if($line2 =~ /^ATOM/){
					$atom   = substr($line2, 13, 3);
					$cid    = substr($line2, 21, 1);
					$resnum = substr($line2, 22, 5); # for insertion code
					$x      = substr($line2, 30, 8);
					$y      = substr($line2, 38, 8);
					$z      = substr($line2, 46, 8);

					if($cid eq "H" && $atom eq "CA " && ( $resnum eq " $base_start " || $resnum eq " 101 " || $resnum eq " 102 " || $resnum eq " 103 " || $resnum eq "104 " || $resnum eq " 105 ") ){
						print DATA "$resnum\t$x\t$y\t$z\t1.00\t0.00\n";
					}

					# Calc distance here	Ca(n+1) - Ca(n-2)
					if($cid eq "H" && $resnum eq " $base_start " && $atom eq "CA "){
						$x_o = substr($line2, 30, 8);
						$y_o = substr($line2, 38, 8);
						$z_o = substr($line2, 46, 8);
					}elsif($cid eq "H" && $resnum eq " 103 " && $atom eq "CA "){
						$x_ne1 = substr($line2, 30, 8);
						$y_ne1 = substr($line2, 38, 8);
						$z_ne1 = substr($line2, 46, 8);

						$pre_dist = ($x_ne1 - $x_o) * ($x_ne1 - $x_o) + ($y_ne1 - $y_o) * ($y_ne1 - $y_o) + ($z_ne1 - $z_o) * ($z_ne1 - $z_o);
						$dist = sqrt($pre_dist);

						print DIST "$dist\n";
					}
				}
			}
			close(PDB);
			close(DATA);
			close(DIST);
		}else{
			print "No PDB...$path/$target\t$num\n";
			#next;
		}

		if(-e "./out_classification/angle/$target\_$num.base"){
			`./calcDihed < ./out_classification/angle/$target\_$num.base >> out_classification/base/$target.base`;
		}else{
			`echo 0 >> out_classification/base/$target.base`
		}
	}


	open(BASE, "out_classification/base/$target.base");
	while(my $line3 = <BASE>){
		$cnt++;

		if($cnt < 10){
			$num = "000$cnt";
		}elsif($cnt < 100){
			$num = "00$cnt";
		}elsif($cnt < 1000){
			$num = "0$cnt";
		}elsif($cnt < 10000){
			$num = "$cnt";
		}

		$line3 =~ s/\n//;
		if($line3 >= -20 && $line3 <= 80){
			$base = "Kink";
			$cnt_k++;
		}elsif($line3 >= 125 || $line3 <= -125){
			$base = "Extend";
			$cnt_e++;
		}else{
			$base = "Other";
			$cnt_o++;
		}

		print ALL "$target\t$num\t$base\t$line3\t";

		if(-e "./out_classification/dist/$target\_$num.dist"){
			open(DIST2, "./out_classification/dist/$target\_$num.dist") or die "check7\t$num\n";
		}else{
			print ALL "0\n";
			next;
		}
		while(my $line4 = <DIST2>){
			$line4 =~ s/\n//;

			print ALL "$line4\n";
		}
		close(DIST2);

	}
	close(BASE);

	print RESULT "$target\t$len_h3\t$cnt_k\t$cnt_e\t$cnt_o\n";
	#exit;

	$cnt_k = 0;
	$cnt_e = 0;
	$cnt_o = 0;
	$cnt   = 0;

	close(ALL);
}
close(LIST);
close(RESULT);

