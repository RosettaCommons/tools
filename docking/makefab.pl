#!/usr/bin/perl

if($#ARGV < 1) {

    print "Usage $0 <pdbfile> <chainids>\n Eg: $0 1AHW.pdb LH\n <L> is light chain id, <H> is heavy chain id\n";
    print " The pdb file name in the fab file is the first four characters of the input pdb file name\n";
    die();
}



%threetoone=("ALA",A,"CYS",C,"ASP",D,"GLU",E,"PHE",F,"GLY",G,"HIS",H,"ILE",I,"LYS",K,"LEU",L,"MET",M,"ASN",N,"PRO",P,"GLN",Q,"ARG",R,"SER",S,"THR",T,"VAL",V,"TRP",W,"TYR",Y);

$l1std=11;
$l2std=7;
$l3std=9;
$h1std=10;
$h2std=16;
$h3std=8;



$pdbfile = shift @ARGV;
$chains = shift @ARGV;

&initialize;
&readpdbfile;
&findcdrs;
&assignfab;
&initialize;

sub initialize{

#get the pdbfilename and the chain ID's 

$firstfour=substr($pdbfile,0,4);
$lightchain=substr($chains,0,1);
$heavychain=substr($chains,1,1);
$nreslight=$ncyslight=0;
$nresheavy=$ncysheavy=0;
$lightseq=$heavyseq=$fablightseq=$fabheavyseq=$l1=$l2=$l3=$h1=$h2=$h3='';

}

sub readpdbfile{

#read the actual pdb file and then identify the light and the heavy chain sequence.

    
    $filename=$pdbfile;

    open(pdbfile,$filename)||die();
    $line=<pdbfile>;
    chop($line);


    while($line ne ''){


	($identifier,$atomno,$atom,$residue,$chain,$residueno,@junk)=split(/ +/,$line) ;

	if($identifier eq "ATOM" and $atom eq "CA") {

	    if($chain eq $lightchain) {

		$lightseq=$lightseq.$threetoone{$residue};
		$nreslight++;
		$fablight[$nreslight]=".";
		$reslight[$nreslight]=$residue;

	    }elsif($chain eq $heavychain) {

		$heavyseq=$heavyseq.$threetoone{$residue} ;
		$nresheavy++ ;
		$fabheavy[$nresheavy]="." ;
		$resheavy[$nresheavy]=$residue ;

	    }

	}

	$line=<pdbfile> ;
	chop($line) ;
    }

    $fablight[1]=$fablight[2]="N" ;

}

sub findcdrs{

#*********L1***************************
$var = $lightseq =~/C[A-Z]{1,20}(WYL|WLQ|WFQ|WYQ|WYH|WVQ|WVR|WWQ|WVK)/;
if($var){
    $temp=$&;
    $lenl1=length ($temp)-4;
    $l1=substr($temp,1,$lenl1);
}
#************************************


#***********L3********************
$var = $lightseq =~/C[A-Z]{1,20}(F|W)G[A-Z]G/;
if($var){
    $temp=$&;
    $lenl3=length ($temp)-5;
    $l3=substr($temp,1,$lenl3);
}
#****************************

#**************H1************
$var = $heavyseq =~/C[A-Z]{1,25}(WIRQ|WVRQ|WVKQ|WFKQ|WIRK|WYQQ|WVKH|WAKQ)/;
if($var){
    $temp=$&if($var);
    $lenh1=length ($temp)-8;
    $h1=substr($temp,4,$lenh1);
}

#******************************



#***********H3****************
$var = $heavyseq =~/C[A-Z]{1,27}(F|W)G[A-Z]G/;
if($var){
    $temp=$&;
    $lenh3=length ($temp)-7;
    $h3=substr($temp,3,$lenh3);
    $h3andstem=substr($temp,0,$lenh3+3);
}
#***************************



$l1start= index($lightseq,$l1);
$l1end=$l1start+$lenl1-1;
$l2start=$l1end+16;
$l2end=$l2start+7-1;
$l3start= index($lightseq,$l3);
$l3end=$l3start+$lenl3-1;
$l2=substr($lightseq,$l2start,7);
$lenl2=7;

$h1start = index($heavyseq,$h1);
$h1end=$h1start+$lenh1-1;
$h3start= index($heavyseq,$h3);
$h3end=$h3start+$lenh3-1;



$h2start=$h1end+15;
$h2end=$h3start-33;
$lenh2=$h2end-$h2start+1;
$h2= substr($heavyseq,$h2start,$lenh2);


$frl1=substr($lightseq,0,$l1start);
$lenfrl1=length($frl1);
$frl2=substr($lightseq,$l1end+1,15);
$lenfrl2=length($frl2);
$frl3=substr($lightseq,$l2end+1,$l3start-$l2end-1);
$lenfrl3=length($frl3);
$frl4=substr($lightseq,$l3end+1);
$lenfrl4=length($frl4);
$frh1=substr($heavyseq,0,$h1start);
$lenfrh1=length($frh1);
$frh2=substr($heavyseq,$h1end+1,$h2start-$h1end-1);
$lenfrh2=length($frh2);
$frh3=substr($heavyseq,$h2end+1,$h3start-$h2end-1);
$lenfrh3=length($frh3);
$frh4=substr($heavyseq,$h3end+1);
$lenfrh4=length($frh4);
$seq1=$frh1.$h1.$frh2.$h2.$frh3.$h3;


}


sub assignfab {



    @lightchaincdr=($l1start+1..$l1end+1,$l2start+1..$l2end+1,$l3start+1..$l3end+1);
    @heavychaincdr=($h1start+1..$h1end+1,$h2start+1..$h2end+1,$h3start+1..$h3end+1);
    @flankcdrlight=($l1start-1,$l1start,$l1end+2,$l1end+3,$l2start-1,$l2start,$l2end+2,$l2end+3,$l3start-1,$l3start,$l3end+2,$l3end+3);
    @flankcdrheavy=($h1start-1,$h1start,$h1end+2,$h1end+3,$h2start-1,$h2start,$h2end+2,$h2end+3,$h3start-1,$h3start,$h3end+2,$h3end+3);


    foreach $position(@lightchaincdr){
	$fablight[$position]="T";
    }

    foreach $position(@heavychaincdr){
	$fabheavy[$position]="T";
    }


    foreach $position(@flankcdrlight){
	$fablight[$position]="N";
    }

    foreach $position(@flankcdrheavy){
	$fabheavy[$position]="N";
    }


    for ($i=1;$i<=$nreslight;$i++){

	$fablightseq=$fablightseq.$fablight[$i];

    }

    for ($i=1;$i<=$nresheavy;$i++){

	$fabheavyseq=$fabheavyseq.$fabheavy[$i];

    }




    $ncharsperline=60;

    $quo = int(($nreslight+1)/$ncharsperline);
    $remain=($nreslight+1)-($ncharsperline*$quo);
    $nlines=$quo;
    $nlines++ if($remain != 0);

    for($i=1;$i<=$nlines;$i++){


	$skip=($i-1)*$ncharsperline;
	$seq=substr($lightseq,$skip,$ncharsperline);
	$fab=substr($fablightseq,$skip,$ncharsperline);
	print "$firstfour $lightchain          $seq\n";
	print "TFN             $fab\n";


    }

    $quo = int(($nresheavy+1)/$ncharsperline);
    $remain=($nresheavy+1)-($ncharsperline*$quo);
    $nlines=$quo;
    $nlines++ if($remain != 0);

    for($i=1;$i<=$nlines;$i++){


	$skip=($i-1)*$ncharsperline;
	$seq=substr($heavyseq,$skip,$ncharsperline);
	$fab=substr($fabheavyseq,$skip,$ncharsperline);
	print "$firstfour $heavychain          $seq\n";
	print "TFN             $fab\n";


    }

}

