#!/usr/bin/perl
if ( @ARGV < 1 ) {

print "Not enough arguments. Exiting...
Usage: $0 <listfile>
Script to re-number PDB files contained in <listfile> according to the Chothia scheme
This script should be run after the script chothia_mapping.pl has been run with the same arguments
The Chothia re-numbered file has the _chothia.pdb extension
<listfile> is a list of PDB file names and chain ID's in the following format
1AHW LH:A
1BQL AB:C etc...
The first column has the PDB ID's (of any length) without the .pdb extension
The second column has the light, heavy and antigen chain ID's separated by a :\n\n" ;
die () ;

}

print ("Creating Chothia numbered files...\n");

$nargs=$#ARGV;
$stdoccupancy=1.00;
$stdbfactor=25.00;

$list=shift @ARGV;
open(LIST,$list)||die();


if($nargs == 2){
	$chothialight = shift @ARGV;
	$chothiaheavy = shift @ARGV;
}


&initialize;

while($pdbfile ne ''){
	print "Renumbering $pdbfile ...$lightchain $heavychain\t$agchains\n";
	&readpdbfile;
	&initialize;
}


sub initialize{
	#get the pdbfilename and the chain ID's 
	$line=<LIST>;
	chop($line);
	@fileandchains=split(/ +/,$line);
	$pdbfile=$fileandchains[0];
	$chains=$fileandchains[1];
	@array=split(/:/,$chains);
	$lightchain=substr($array[0],0,1);
	$heavychain=substr($array[0],1,1);

	$agchains=substr($array[1], 0);
	$agchains =~ s/\n//;
}

sub readpdbfile{
	$filename=$pdbfile.".pdb";
#    $filename = $pdbfile.".ent" if !(-e $filename);
	open(PDBFILE,$filename)||die();

	$filename1 = $pdbfile."\_chothia.light";
	if(-e $filename1){
		open(CHOTHIALIGHT,$filename1)||die();
	}else{
		open(CHOTHIALIGHT,$chothialight)||die();
	}


	$filename2 = $pdbfile."\_chothia.heavy";
	if(-e $filename2){
		open(CHOTHIAHEAVY,$filename2)||die();
	}else{
		open(CHOTHIAHEAVY,$chothiaheavy)||die();
	}

	&readchothialight;
	&readchothiaheavy;

	$line=<PDBFILE>;
	chop($line);

	$var1=$var2=0;
	%newlightnumber=%newheavynumber=();

	while($line ne ''){
		$identifier = substr($line,0,6);
		$atomno = substr($line,6,5);
		$atom = substr($line,12,4);
		$alt_loc = substr($line,16,1);
		$residue = substr($line,17,3);
		$chain = substr($line,21,1);
		$residueno = substr($line,22,4);
		$insert_code=substr($line,26,1);
		$x = substr($line,30,8);
		$y = substr($line,38,8);
		$z =  substr($line,46,8);
		$bfactor =  &removespaces (substr($line, 60,6)) if ( $identifier =~ /ATOM/ ) ;

		if($identifier =~ "ATOM" and $chain eq $lightchain and $atom =~ "CA"){
			unless(  ($old_residueno eq $residueno) and ($old_alt_loc ne $alt_loc) and ($old_insert_code eq $insert_code)   ){
				$old_alt_loc =$alt_loc;
				$old_residueno=$residueno;

				$var1++;
				$argument=$residueno.$insert_code;
				$newlightnumber{$argument}=$chothialight[$var1];

				$old_alt_loc =$alt_loc;
				$old_residueno=$residueno;
				$old_insert_code=$insert_code;
			}
		}elsif($identifier =~ "ATOM" and $chain eq $heavychain and $atom =~ "CA"){
			unless(   ($old_residueno eq $residueno)   and  ($old_alt_loc ne $alt_loc)  and ($old_insert_code eq $insert_code)   ){
				$old_alt_loc =$alt_loc;
				$old_residueno=$residueno;
				$old_insert_code=$insert_code; 
				$var2++;
				$argument=$residueno.$insert_code;
				$newheavynumber{$argument}=$chothiaheavy[$var2];
			}
		}
		$line=<PDBFILE>;
		chop($line);
	}

	close(PDBFILE);
	open(PDBFILE,$filename)||die();

	$pdbfile =~ s/pdb_tmp\///;

	open(CHOTHIAPDB,">pdb$pdbfile\_chothia.pdb");

	$line=<PDBFILE>;
	chop($line);

	while($line ne ''){
			$identifier = substr($line,0,6);
			$atomno = substr($line,6,5);
			$atom = substr($line,12,4);
			$alt_loc = substr($line,16,1);
			$residue = substr($line,17,3);
			$chain = substr($line,21,1);
			$residueno = substr($line,22,4);
			$insert_code=substr($line,26,1);
			$x = substr($line,30,8);
			$y = substr($line,38,8);
			$z =  substr($line,46,8);
			$occupancy = substr($line,54,6);
			$bfactor =  &removespaces (substr($line, 60,6)) if ( $identifier =~ /ATOM/) ;
			$var=0;

			if($identifier =~ "ATOM" and $chain eq $lightchain){
			$argument=$residueno.$insert_code;
			$residueno_new=$newlightnumber{$argument};

			$insert_code='';
			$var = $residueno_new =~ /[A-Z]/;
			if($var){
					$insert_code=$&;
					$residueno_new=substr($residueno_new,0,length($residueno_new)-1);
			}
			printf CHOTHIAPDB "%6s%5s %-4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f              \n","ATOM  ",$atomno,$atom,$alt_loc,$residue,$chain,$residueno_new,$insert_code,$x,$y,$z,$stdoccupancy,$bfactor if($newlightnumber{$argument} ne '');
			}elsif($identifier =~ "ATOM" and $chain eq $heavychain){
				$argument=$residueno.$insert_code;
				$residueno_new=$newheavynumber{$argument};
				#$residueno_new=$newheavynumber{$residueno};

				$insert_code='';
				$var = $residueno_new =~ /[A-Z]/;

				if($var){
					$insert_code=$&;
					$residueno_new=substr($residueno_new,0,length($residueno_new)-1);
				}

				printf CHOTHIAPDB "%6s%5s %-4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f              \n","ATOM  ",$atomno,$atom,$alt_loc,$residue,$chain,$residueno_new,$insert_code,$x,$y,$z,$stdoccupancy,$bfactor if($newheavynumber{$argument} ne '');
			}elsif($identifier =~ "ATOM" and $agchains =~ /$chain/){
				print CHOTHIAPDB "$line\n";
			}elsif($identifier =~ "HEADER" || $identifier =~ "TITLE" || $line =~ /PH          / || $identifier =~ "COMPND" || $identifier =~ "SOURCE" || $line =~ /RESOLUTION\.   / || $line =~ /R VALUE            \(WORKING SET\)/|| $line =~ /   FREE R VALUE    /){
				print CHOTHIAPDB "$line\n";
			}elsif($identifier =~ "HETATM" && $line =~ /HOH/){
				#print CHOTHIAPDB "$line\n";
			}

		$line=<PDBFILE>;
		chop($line);
	}
}

sub readchothialight{
	@chothialight=[];
	$var1=0;
	$temp=<CHOTHIALIGHT>;
	chop($temp);

	while($temp ne ''){
		$var1++;
		$chothialight[$var1]=$temp;
		$temp=<CHOTHIALIGHT>;
		chop($temp);
	}
}

sub readchothiaheavy{
	@chothiaheavy=[];
	$var2=0;
	$temp=<CHOTHIAHEAVY>;
	chop($temp);

	while($temp ne ''){
		$var2++;
		$chothiaheavy[$var2]=$temp;
		$temp=<CHOTHIAHEAVY>;
		chop($temp);
	}
}

sub removespaces {
	$argument = $_[0] ;
	$argument =~ s/\s+$// ;
	$argument =~ s/^\s+// ;

	return $argument ;
}


