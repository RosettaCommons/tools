#!/usr/bin/perl -w

use strict;

if ( @ARGV < 1 ) {
	print "Not enough arguments. Exiting...
	Usage: $0 <listfile>
	Script to map the given PDB sequence to the Chothia numbering
	The script outputs the .light and .heavy files that contain the mapping
	<listfile> is a list of PDB file names and chain ID's in the following format
	1AHW LH:A
	1BQL AB:C etc...
	The first column has the PDB ID's (of any length) without the .pdb extension
	The second column has the light, heavy and antigen chain ID's separated by a :\nMake sure that if the CDR L2's length is not seven (maybe eleven) then you will have to use the chothia_mapping.pl.L211 instead and make changes accordingly\n\n" ;

	die () ;
}

my %threetoone=("ALA","A","CYS","C","ASP","D","GLU","E","PHE","F","GLY","G","HIS","H","ILE","I","LYS","K","LEU","L","MET","M",
                "ASN","N","PRO","P","GLN","Q","ARG","R","SER","S","THR","T","VAL","V","TRP","W","TYR","Y");
# dead code
#my $l1std=11;
#my $l2std=7;
#my $l3std=9;
#my $h1std=10;
#my $h2std=16;
#my $h3std=8;

my ($lenl1,$lenl2,$lenl3,$lenh1,$lenh2,$lenh3);
my ($lenfrl1,$lenfrl2,$lenfrl3,$lenfrl4,$lenfrh1,$lenfrh2,$lenfrh3,$lenfrh4);
my ($frl1,$frl2,$frl3,$frl4,$frh1,$frh2,$frh3,$frh4);

my ($pdbfile,$filename);
my ($lightchain, $heavychain);
my ($nreslight,$ncyslight);
my ($nresheavy,$ncysheavy);
my ($lightseq,$heavyseq);
my ($l1,$l2,$l3,$h1,$h2,$h3);

my (@newnumberfrl1, @newnumberfrl2, @newnumberfrl3, $newnumberfrl4, @newnumberl1, @newnumberl2, @newnumberl3);
my (@newnumberfrh1, @newnumberfrh2, @newnumberfrh3, $newnumberfrh4, @newnumberh1, @newnumberh2, @newnumberh3);

my ($heavyseq_first, $heavyseq_second);
my ($lightseq_first, $lightseq_second);

my $len;


my $list=shift @ARGV;
open(LIST,"<$list")||die "E: Could not open '$list'. $!";

&initialize;

while($pdbfile ne ''){
	print "Mapping $pdbfile ...\n";

	&readpdbfile;
	&findcdrs;
	&assignnumbering;
	&checknumbering;
	&renumbercdrs;
	&initialize;
}

sub initialize{
	#get the pdbfilename and the chain ID's 
	my $line=<LIST>;
	if (defined($line)) {
		chomp($line);
		my @fileandchains=split(/ +/,$line);
		$pdbfile=$fileandchains[0];
		my $chains=$fileandchains[1];
		my @array=split(/:/,$chains);
		$lightchain=substr($array[0],0,1);
		$heavychain=substr($array[0],1,1);
		$nreslight=$ncyslight=0;
		$nresheavy=$ncysheavy=0;
		$lightseq=$heavyseq=$l1=$l2=$l3=$h1=$h2=$h3='';
	} else {
		# terminates the execution
		$pdbfile="";
	}
}

sub readpdbfile {
	#read the actual pdb file and then identify the light and the heavy chain sequence.
	$filename=$pdbfile.".pdb";

	#$filename = $pdbfile.".ent" if !(-e $filename);

	unless( -e $filename){
		print "No PDB $filename\n";
		exit;
	}


	open(PDBFILE,"<$filename") or die("Could not open pdfile '$filename'. $!");
	open(CHOTHIALIGHT,">$pdbfile\_chothia.light") or die ("Could not create cothia-numbered light chain at '$pdbfile\_chothia.light'. $!");
	open(CHOTHIAHEAVY,">$pdbfile\_chothia.heavy") or die ("Could not create cothia-numbered heavy chain at '$pdbfile\_chothia.heavy'. $!");
	my $line=<PDBFILE>;
	chop($line);

	while($line ne ''){
		my ($identifier,$atomno,$atom,$residue,$chain,$residueno,@junk)=split(/ +/,$line);

		$identifier = substr($line,0,6);
		$atomno = substr($line,6,5);
		$atom = substr($line,12,4);
		my $alt_loc = substr($line,16,1);
		$residue = substr($line,17,3);
		$chain = substr($line,21,1);
		$residueno = substr($line,22,4);
		my $insert_code=substr($line,26,1);
		my $x = substr($line,30,8);
		my $y = substr($line,38,8);
		my $z = substr($line,46,8);

		my ($old_residueno,$old_insert_code,$old_alt_loc)="";
		if($identifier =~ "ATOM" and $atom =~ "CA"){
			if($chain eq $lightchain){
				$lightseq=$lightseq.$threetoone{$residue}unless(   ($old_residueno eq $residueno) and ($old_alt_loc ne $alt_loc)  and ($old_insert_code eq $insert_code)  );
				$old_residueno=$residueno;
				$old_alt_loc=$alt_loc;
				$old_insert_code=$insert_code
			}elsif($chain eq $heavychain){
				$heavyseq=$heavyseq.$threetoone{$residue} unless(($old_residueno eq $residueno) and ($old_alt_loc ne $alt_loc)  and ($old_insert_code eq $insert_code) );
				$old_residueno=$residueno;
				$old_alt_loc=$alt_loc;
				$old_insert_code=$insert_code
			}
		}

		$line=<PDBFILE>;
		if (defined($line)) {
			chomp($line);
		} else {
			$line="";
		}
	}

	### Split a sequence into two parts ###
	if(length($lightseq) > 120){
		$lightseq_first  = substr($lightseq, 0, 60);
		$lightseq_second = substr($lightseq, 50, 70);
	}else{
		$lightseq_first  = substr($lightseq, 0, 60);
		$lightseq_second = substr($lightseq, 50);
	}

	if(length($heavyseq) > 140){
		$heavyseq_first  = substr($heavyseq, 0, 60);
		$heavyseq_second = substr($heavyseq, 50, 90);
	}else{
		$heavyseq_first  = substr($heavyseq, 0, 60);
		$heavyseq_second = substr($heavyseq, 50);
	}
}

sub findcdrs{
	#*********L1***************************
	# C[A-Z]{1,17}(WYL|WLQ|WFQ|WYQ|WYH|WVQ|WVR|WWQ|WVK|WYR|WLL|WFL|WVF|WIQ|WYR|WNQ|WHL|WHQ|WYM|WYY)
	my $var;
	$var = $lightseq_first =~/C[A-Z]{1,17}(WYL|WLQ|WFQ|WYQ|WYH|WVQ|WVR|WWQ|WVK|WYR|WLL|WFL|WVF|WIQ|WYR|WNQ|WHL|WYM|WYY)/;
	if($var){
		my $temp=$&;
		$lenl1=length ($temp)-4;
		$l1=substr($temp,1,$lenl1);
	}
	#************************************

	#***********L3********************
	# C[A-Z]{1,15}(F|V|S)G[A-Z](G|Y)
	$var = $lightseq_second =~/C[A-Z]{1,15}(L|F|V|S)G[A-Z](G|Y)/;
	if($var){
		my $temp=$&;
		$lenl3=length ($temp)-5;
		$l3=substr($temp,1,$lenl3);

	}
	if($pdbfile eq "1xiw"){
		$lenl3 = 9;
		$l3 = "QQGNTLPWT";
	}
	#****************************

	#**************H1************
	# C[A-Z]{1,16}(W)(I|V|F|Y|A|M|L|N|G)(R|K|Q|V|N|C|G)(Q|K|H|E|L|R)
	$var = $heavyseq_first =~/C[A-Z]{1,16}(W)(I|V|F|Y|A|M|L|N|G)(R|K|Q|V|N|C)(Q|K|H|E|L|R)/; 
	if($var){
		my $temp=$&;
		$lenh1=length ($temp)-8;
		$h1=substr($temp,4,$lenh1);
	}
	#******************************

	$var = 0;

	#***********H3****************
	if($pdbfile eq "1tqb"){
		$var = $heavyseq_second =~/FTRGTDYWGQG/; # for 1ghf:"WAQG", for 3se8:"WCQG", for 3mug,3u1s,3u2s:more than 27
		if($var){
			my $temp=$&;
			$lenh3=length ($temp)-7;
			$h3=substr($temp,3,$lenh3);
			$h1 = "GYTFTNYGMN";
			$lenh1 = 10;
		}
	}else{
		# C[A-Z]{1,33}(W)(G|A|C)[A-Z](S|G|R)
		$var = $heavyseq_second =~/C[A-Z]{1,33}(W)(G|A|C)[A-Z](S|Q|G|R)/; # for 1ghf:"WAQG", for 3se8:"WCQG", for 3mug,3u1s,3u2s:more than 27
		if($var){
			my $temp=$&;
			$lenh3=length ($temp)-7;
			$h3=substr($temp,3,$lenh3);
		}
	}

	#print "$pdbfile\t$var\t$heavyseq_second\t$temp\n";

	#***************************

	my $l1start= index($lightseq,$l1);
	my $l1end=$l1start+$lenl1-1;
	my $l2start=$l1end+16;
	my $l2end=$l2start+7-1;
	#L211
	#$l2end=$l2start+11-1;
	my $l3start= index($lightseq,$l3);
	my $l3end=$l3start+$lenl3-1;
	my $l2=substr($lightseq,$l2start,7);
	$lenl2=7;
	#L211
	#$lenl2=11;

	my $h1start = index($heavyseq,$h1);
	my $h1end=$h1start+$lenh1-1;
	my $h3start= index($heavyseq,$h3);
	my $h3end=$h3start+$lenh3-1;

	my $h2start=$h1end+15;
	my $h2end=$h3start-33;

	if($pdbfile eq "1MFE"){
		$h2start = $h1end + 15 - 1;
	}elsif($pdbfile eq "1MRD" || $pdbfile eq "1MRF"){
		$h2end = $h3start - 33 + 3;
	}elsif($pdbfile eq "2X7L"){
		$h2end = $h3start - 33 + 2;
	}elsif($pdbfile eq "3MNW"){
		$h2start = $h1end + 15 - 2;
	}


	$lenh2=$h2end-$h2start+1;
	$h2= substr($heavyseq,$h2start,$lenh2);

	$frl1=substr($lightseq,0,$l1start);
	$lenfrl1=length($frl1);
	$frl2=substr($lightseq,$l1end+1,15);
	$lenfrl2=length($frl2);
	$frl3=substr($lightseq,$l2end+1,$l3start-$l2end-1);
	$lenfrl3=length($frl3);
	$frl4=substr($lightseq,$l3end+1,12);
	$lenfrl4=length($frl4);
	$frh1=substr($heavyseq,0,$h1start);
	$lenfrh1=length($frh1);
	$frh2=substr($heavyseq,$h1end+1,$h2start-$h1end-1);
	$lenfrh2=length($frh2);
	$frh3=substr($heavyseq,$h2end+1,$h3start-$h2end-1);
	$lenfrh3=length($frh3);
	$frh4=substr($heavyseq,$h3end+1,10);
	$lenfrh4=length($frh4);
	my $seq1=$frh1.$h1.$frh2.$h2.$frh3.$h3; # not used again

	print "$filename\t$frl1 - \"$l1\" - $frl2 - \"$l2\" - $frl3 - \"$l3\" - $frl4\n";
	print "$filename\t$frh1 - \"$h1\" - $frh2 - \"$h2\" - $frh3 - \"$h3\" - $frh4\n";
}


sub renumbercdrs{
	my @string;
	$string[1]=$newnumberfrl1[$lenfrl1];
	$string[2]=$newnumberl1[$lenl1];
	$string[3]=$newnumberfrl2[$lenfrl2];
	$string[4]=$newnumberl2[$lenl2];
	$string[5]=$newnumberfrl3[$lenfrl3];
	$string[6]=$newnumberl3[$lenl3];
	$string[7]=$newnumberfrl4;

	for(my $i=1;$i<=7;$i++){
		my @array=split(/,/,$string[$i]);
		my $nelements=@array;

		for(my $j=0;$j <$nelements;$j++){
			print CHOTHIALIGHT "$array[$j]\n";
		}

	}

	$string[1]=$newnumberfrh1[$lenfrh1];
	$string[2]=$newnumberh1[$lenh1];
	$string[3]=$newnumberfrh2[$lenfrh2];
	$string[4]=$newnumberh2[$lenh2];
	$string[5]=$newnumberfrh3[$lenfrh3];
	$string[6]=$newnumberh3[$lenh3];
	$string[7]=$newnumberfrh4;

	for(my $i=1;$i<=7;$i++){
		my @array=split(/,/,$string[$i]);
		my $nelements=@array;

		for(my $j=0;$j <$nelements;$j++){
			print CHOTHIAHEAVY "$array[$j]\n";
		}
	}
}

sub assignnumbering{
	$newnumberfrl1[20]="3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23";
	$newnumberfrl1[21]="2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23";
	$newnumberfrl1[22]="1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23";
	$newnumberfrl1[23]="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23";
	$newnumberfrl1[24]="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23";
	$newnumberfrl1[25]="-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23";
	$newnumberfrl1[26]="-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23";

	if($frl1 =~ /[LMVI][A-Z][QE][A-Z]{9}G[A-Z]{4}[LVIMF][STN]C/){
		$newnumberfrl1[20]="4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23";
		$newnumberfrl1[21]="3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23";
		$newnumberfrl1[22]="2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23";
		$newnumberfrl1[23]="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23";
		$newnumberfrl1[24]="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23";
	}elsif ($frl1 =~ /[LMVI][A-Z][QE][A-Z]{8}G[A-Z]{4}[LVIMF][STN]C/){
		$newnumberfrl1[20]="3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23";
		$newnumberfrl1[21]="2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23";
		$newnumberfrl1[22]="1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23";
		$newnumberfrl1[23]="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23";
		$newnumberfrl1[24]="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23";
	}

	$newnumberl1[8] ="24,25,26,27,28,29,30,34";
	$newnumberl1[9] ="24,25,26,27,28,29,30,33,34";
	$newnumberl1[10]="24,25,26,27,28,29,30,32,33,34";
	$newnumberl1[11]="24,25,26,27,28,29,30,31,32,33,34";
	$newnumberl1[12]="24,25,26,27,28,29,30,30A,31,32,33,34";
	$newnumberl1[13]="24,25,26,27,28,29,30,30A,30B,31,32,33,34";
	$newnumberl1[14]="24,25,26,27,28,29,30,30A,30B,30C,31,32,33,34";
	$newnumberl1[15]="24,25,26,27,28,29,30,30A,30B,30C,30D,31,32,33,34";
	$newnumberl1[16]="24,25,26,27,28,29,30,30A,30B,30C,30D,30E,31,32,33,34";
	$newnumberl1[17]="24,25,26,27,28,29,30,30A,30B,30C,30D,30E,30F,31,32,33,34";
	
	$newnumberfrl2[15]="35,36,37,38,39,40,41,42,43,44,45,46,47,48,49";
	
	$newnumberl2[7]="50,51,52,53,54,55,56";
	$newnumberl2[8]="50,51,52,53,54,54A,55,56";
	$newnumberl2[9]="50,51,52,53,54,54A,54B,55,56";
	$newnumberl2[10]="50,51,52,53,54,54A,54B,54C,55,56";
	$newnumberl2[11]="50,51,52,53,54,54A,54B,54C,54D,55,56";
	
	# Add length 35/36
	$newnumberfrl3[32]="57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88";
	$newnumberfrl3[33]="57,58,59,60,61,62,63,64,65,66,66A,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88";
	$newnumberfrl3[34]="57,58,59,60,61,62,63,64,65,66,66A,66B,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88";
	$newnumberfrl3[35]="57,58,59,60,61,62,63,64,65,66,66A,66B,66C,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88";
	$newnumberfrl3[36]="57,58,59,60,61,62,63,64,65,66,66A,66B,66C,66D,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88";

	# Add length 5,6,7	
	$newnumberl3[5] ="89,90,91,92,97"; ## L89 - L97, L98
	$newnumberl3[6] ="89,90,91,92,93,97";
	$newnumberl3[7] ="89,90,91,92,93,94,97";
	$newnumberl3[8] ="89,90,91,92,93,94,95,97";
	$newnumberl3[9] ="89,90,91,92,93,94,95,96,97";
	$newnumberl3[10]="89,90,91,92,93,94,95,95A,96,97";
	$newnumberl3[11]="89,90,91,92,93,94,95,95A,95B,96,97";
	$newnumberl3[12]="89,90,91,92,93,94,95,95A,95B,95C,96,97";
	$newnumberl3[13]="89,90,91,92,93,94,95,95A,95B,95C,95D,96,97";
	$newnumberl3[14]="89,90,91,92,93,94,95,95A,95B,95C,95D,95E,96,97";
	$newnumberl3[15]="89,90,91,92,93,94,95,95A,95B,95C,95D,95E,95F,96,97";
	
	$newnumberfrl4="98,99,100,101,102,103,104,105,106,107,108,109";
	
	$newnumberfrh1[21]="5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25";
	$newnumberfrh1[22]="4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25";
	$newnumberfrh1[23]="3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25";
	$newnumberfrh1[24]="2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25";
	$newnumberfrh1[25]="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25";
	$newnumberfrh1[26]="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25";
	
	# Add length 6
	$newnumberh1[6]="26,27,32,33,34,35";
	$newnumberh1[7]="26,27,28,32,33,34,35";
	$newnumberh1[8]="26,27,28,29,32,33,34,35";
	$newnumberh1[9]="26,27,28,29,30,32,33,34,35";
	$newnumberh1[10]="26,27,28,29,30,31,32,33,34,35";
	$newnumberh1[11]="26,27,28,29,30,31,31A,32,33,34,35";
	$newnumberh1[12]="26,27,28,29,30,31,31A,31B,32,33,34,35";
	$newnumberh1[13]="26,27,28,29,30,31,31A,31B,31C,32,33,34,35";
	
	$newnumberfrh2[14]="36,37,38,39,40,41,42,43,44,45,46,47,48,49";
	
	$newnumberh2[10]="50,51,52,59,60,61,62,63,64,65";
	$newnumberh2[11]="50,51,52,58,59,60,61,62,63,64,65";
	$newnumberh2[12]="50,51,52,57,58,59,60,61,62,63,64,65";
	$newnumberh2[13]="50,51,52,56,57,58,59,60,61,62,63,64,65";
	$newnumberh2[14]="50,51,52,55,56,57,58,59,60,61,62,63,64,65";
	$newnumberh2[15]="50,51,52,54,55,56,57,58,59,60,61,62,63,64,65";
	$newnumberh2[16]="50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65";
	$newnumberh2[17]="50,51,52,52A,53,54,55,56,57,58,59,60,61,62,63,64,65";
	$newnumberh2[18]="50,51,52,52A,52B,53,54,55,56,57,58,59,60,61,62,63,64,65";
	$newnumberh2[19]="50,51,52,52A,52B,52C,53,54,55,56,57,58,59,60,61,62,63,64,65";
	$newnumberh2[20]="50,51,52,52A,52B,52C,52D,53,54,55,56,57,58,59,60,61,62,63,64,65";
	$newnumberh2[21]="50,51,52,52A,52B,52C,52D,52E,53,54,55,56,57,58,59,60,61,62,63,64,65";
	$newnumberh2[22]="50,51,52,52A,52B,52C,52D,52E,52F,53,54,55,56,57,58,59,60,61,62,63,64,65";
	
	# Add length 29/30/31
	$newnumberfrh3[29]="66,67,68,69,70,71,72,76,77,78,79,80,81,82,82A,82B,82C,83,84,85,86,87,88,89,90,91,92,93,94";
	$newnumberfrh3[30]="66,67,68,69,70,71,72,73,76,77,78,79,80,81,82,82A,82B,82C,83,84,85,86,87,88,89,90,91,92,93,94";
	$newnumberfrh3[31]="66,67,68,69,70,71,72,73,74,76,77,78,79,80,81,82,82A,82B,82C,83,84,85,86,87,88,89,90,91,92,93,94";
	$newnumberfrh3[32]="66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,82A,82B,82C,83,84,85,86,87,88,89,90,91,92,93,94";
	
	$newnumberh3[3]="95,101,102";
	$newnumberh3[4]="95,96,101,102";
	$newnumberh3[5]="95,96,97,101,102";
	$newnumberh3[6]="95,96,97,98,101,102";
	$newnumberh3[7]="95,96,97,98,99,101,102";
	$newnumberh3[8]="95,96,97,98,99,100,101,102";
	$newnumberh3[9]="95,96,97,98,99,100,100A,101,102";
	$newnumberh3[10]="95,96,97,98,99,100,100A,100B,101,102";
	$newnumberh3[11]="95,96,97,98,99,100,100A,100B,100C,101,102";
	$newnumberh3[12]="95,96,97,98,99,100,100A,100B,100C,100D,101,102";
	$newnumberh3[13]="95,96,97,98,99,100,100A,100B,100C,100D,100E,101,102";
	$newnumberh3[14]="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,101,102";
	$newnumberh3[15]="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,101,102";
	$newnumberh3[16]="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,101,102";
	$newnumberh3[17]="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,101,102";
	$newnumberh3[18]="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,101,102";
	$newnumberh3[19]="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,101,102";
	$newnumberh3[20]="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,101,102";
	$newnumberh3[21]="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,101,102";
	$newnumberh3[22]="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,101,102";
	$newnumberh3[23]="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,101,102";
	$newnumberh3[24]="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,101,102";
	$newnumberh3[25]="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,101,102";
	$newnumberh3[26]="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,101,102";
	$newnumberh3[27]="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,101,102";
	$newnumberfrh4="103,104,105,106,107,108,109,110,111,112";
}


sub checknumbering{
	my $error="NOK";
	my ($error1,$error2,$error3,$error4,$error5)=(0,0,0,0,0);

	my $val1=substr($frh1,-22,1);
	my $val2=substr($frh1,-20,1);
	my $val3=substr($frh1,-11,1);
	my $val4=substr($frh1,-6,1);
	my $val5=substr($frh1,-5,1);

	my $val6=substr($frh1,-19,1);
	my $val7=substr($frh1,-18,1);
	my $val8=substr($frh1,-17,1);

	my $h18 = substr($frh1,-8,1);
	my $string1="$val2$val7$val8";
	my $string2="$val2$val6$val7";

	my $type = "unknown";
	$type = "type1" if($string1 eq "EGP");
	$type = "type2" if($string1 eq "EGG");
	$type = "type3" if ($string2 =~ /Q[A-O,Q-Z]G/);
	$type = "type4" if ($string2 eq "QPG");

	$error1=1 if ($val1 =~ /[LVIM]/);
	$error2=1 if ($val2 =~ /[E]/);
	$error3=1 if ($val3 =~ /[GS]/);
	$error4=1 if ($val4 =~ /[LVIMF]/);
	$error5=1 if ($val5 =~ /[ST]/);

	$error ="OK" if($error1 and $error2 and $error3 and $error4 and $error5);
	#print "$pdbfile\t$val1\t$val2\t$val3\t$val4\t$val5\t$error\t$lenfrh1\t$frh1\t$lenfrh3\t$frh3\n";
	#print "DK:$pdbfile\t$lenfrh3\t$frh3\n";

	my $krithi=0; # not used anywhere else
	if($krithi){
		$val1=substr($frl1,-20,1);
		$val2=substr($frl1,-18,1);
		$val3=substr($frl1,-8,1);
		$val4=substr($frl1,-3,1);
		$val5=substr($frl1,-2,1);

		$error1=1 if ($val1 =~ /[LVIM]/);
		$error2=1 if ($val2 =~ /[QE]/);
		$error3=1 if ($val3 =~ /[GS]/);
		$error4=1 if ($val4 =~ /[LVIMF]/);
		$error5=1 if ($val5 =~ /[STN]/);
	}

	$len=length($frl1);
}
