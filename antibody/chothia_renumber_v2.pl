#!/usr/bin/perl -w

use strict;

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

my $nargs=$#ARGV;
my $stdoccupancy=1.00;
my $stdbfactor=25.00;

my $pdbfile;
my ($lightchain, $heavychain, $agchains);

my $list=shift @ARGV;
open(LIST,$list)||die("Could not open file '$list': $!\n");


my (@chothialight, @chothiaheavy);
my ($chothialight, $chothiaheavy); # filenames
if ($nargs == 2) {
	$chothialight = shift @ARGV;
	$chothiaheavy = shift @ARGV;
}

my %newlightnumber;
my %newheavynumber;

my $debug=0;

&initialize;

while($pdbfile ne ''){
	print "Renumbering $pdbfile ...$lightchain $heavychain\t$agchains\n";
	&readpdbfile;
	&initialize;
}


sub initialize{
	#get the pdbfilename and the chain ID's 
	my $line=<LIST>;
	unless(defined($line)) {
		$pdbfile=""; # ends loop
		$line="";
	} else {
		chomp($line);
		my @fileandchains=split(/ +/,$line);
		$pdbfile=$fileandchains[0];
		my $chains=$fileandchains[1];
		my @array=split(/:/,$chains);
		$lightchain=substr($array[0],0,1);
		$heavychain=substr($array[0],1,1);

		$agchains=substr($array[1], 0);
		$agchains =~ s/\n//;
	}
}

sub readpdbfile{
	my $filename=$pdbfile.".pdb";
#    $filename = $pdbfile.".ent" if !(-e $filename);
	open(PDBFILE,"<$filename") or die("Could not open file '$filename': $!\n");

	my $filename1 = $pdbfile."\_chothia.light";
	if(-e $filename1){
		open(CHOTHIALIGHT,"<$filename1") or die("Could not open file '$filename1': $!\n");
	}else{
		open(CHOTHIALIGHT,"<$chothialight") or die("Could not open file '$chothialight': $!\n");
	}


	my $filename2 = $pdbfile."\_chothia.heavy";
	if(-e $filename2){
		open(CHOTHIAHEAVY,"<$filename2") or die("Could not open file '$filename2: $!\n");
	}else{
		open(CHOTHIAHEAVY,"<$chothiaheavy") or die("Could not open file '$chothiaheavy: $!\n");
	}

	&readchothialight;
	&readchothiaheavy;

	my ($var1,$var2)=(0,0);
	my ($old_residueno,$old_alt_loc,$old_insert_code)=("","","");
	while(<PDBFILE>){
		chomp;
		my $line=$_;
		#print "$line\n";
		my $identifier = substr($line,0,6);
		my $atomno = substr($line,6,5);
		my $atom = substr($line,12,4);
		my $alt_loc = substr($line,16,1);
		my $residue = substr($line,17,3);
		my $chain = substr($line,21,1);
		my $residueno = &removespaces(substr($line,22,4));
		my $insert_code=&removespaces(substr($line,26,1));
		my $x = substr($line,30,8);
		my $y = substr($line,38,8);
		my $z = substr($line,46,8);
		my $bfactor = &removespaces (substr($line, 60,6)) if ( $identifier =~ /ATOM/ ) ;

		if($identifier =~ "ATOM" and $chain eq $lightchain and $atom =~ "CA"){
			unless(  ($old_residueno eq $residueno) and ($old_alt_loc ne $alt_loc) and ($old_insert_code eq $insert_code)   ){
				$old_alt_loc =$alt_loc;
				$old_residueno=$residueno;

				$var1++;
				next unless defined ($chothialight[$var1]);
				my $argument=$residueno.$insert_code;
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
				next unless defined ($chothiaheavy[$var2]);
				my $argument=$residueno.$insert_code;
				$newheavynumber{$argument}=$chothiaheavy[$var2];
			}
		}
	}
	close(PDBFILE);

	open(PDBFILE,"<$filename") or die("Could not open file '$filename': $!\n");
	open(CHOTHIAPDB,">pdb$pdbfile\_chothia.pdb") or die("Could not open for writing 'pdb$pdbfile\_chothia.pdb': $!\n");

	while(<PDBFILE>){
		chomp;
		my $line=$_;
		my $identifier = substr($line,0,6);
		my $atomno = substr($line,6,5);
		my $atom = substr($line,12,4);
		my $alt_loc = substr($line,16,1);
		my $residue = substr($line,17,3);
		my $chain = substr($line,21,1);
		my $residueno = &removespaces(substr($line,22,4));
		my $insert_code=&removespaces(substr($line,26,1));
		my $x = substr($line,30,8);
		my $y = substr($line,38,8);
		my $z = substr($line,46,8);
		my $occupancy = substr($line,54,6);
		my $bfactor =  &removespaces (substr($line, 60,6)) if ( $identifier =~ /ATOM/) ;
		my $var=0;

		if($identifier =~ "ATOM" and $chain eq $lightchain){
			my $argument=$residueno.$insert_code;
			next unless exists($newlightnumber{$argument}) and $newlightnumber{$argument} ne '';
			unless (1 or exists($newlightnumber{$argument})) {
				print STDERR "Could not find argument '$argument' among keys of '\$newlightnumber'.\n";
				print STDERR "These are the keys: \n";
				print STDERR "'".join("','",sort keys %newlightnumber)."'\n";
				exit 1;
			}
			my $residueno_new=$newlightnumber{$argument};

			$insert_code='';
			$var = $residueno_new =~ /[A-Z]/;
			if($var){
				$insert_code=$&;
				$residueno_new=substr($residueno_new,0,length($residueno_new)-1);
			}
			printf CHOTHIAPDB "%6s%5s %-4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f              \n",
			                   "ATOM  ",
			                      $atomno,
			                           $atom,
			                              $alt_loc,
			                                 $residue,
			                                     $chain,
			                                        $residueno_new,
			                                           $insert_code,
			                                                  $x,$y,$z,   $stdoccupancy,$bfactor;
		}elsif($identifier =~ "ATOM" and $chain eq $heavychain){
			my $argument=$residueno.$insert_code;
			print STDERR "argument='$argument'\n" if $debug;
			next unless exists($newheavynumber{$argument}) and $newheavynumber{$argument} ne '';
			unless(1 or exists($newheavynumber{$argument})) {
				print STDERR "Could not find argument '$argument' among keys of '\$newheavynumber'.\n";
				print STDERR "These are the keys: \n";
				print STDERR "'".join("','",sort keys %newheavynumber)."'\n";
				exit 1;
			}
			my $residueno_new=$newheavynumber{$argument};
			#$residueno_new=$newheavynumber{$residueno};

			$insert_code='';
			my $var = $residueno_new =~ /[A-Z]/;

			if($var){
				$insert_code=$&;
				$residueno_new=substr($residueno_new,0,length($residueno_new)-1);
			}

			printf CHOTHIAPDB "%6s%5s %-4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f              \n",
			                   "ATOM  ",
			                      $atomno,
			                          $atom,
			                              $alt_loc,
			                                 $residue,
			                                     $chain,
			                                        $residueno_new,
			                                            $insert_code,
			                                                  $x,$y,$z,$stdoccupancy,$bfactor;
		}elsif($identifier =~ "ATOM" and $agchains =~ /$chain/){
				print CHOTHIAPDB "$line\n";
		}elsif($identifier =~ "HEADER" || $identifier =~ "TITLE" || $line =~ /PH          / || $identifier =~ "COMPND" || $identifier =~ "SOURCE" || $line =~ /RESOLUTION\.   / || $line =~ /R VALUE            \(WORKING SET\)/|| $line =~ /   FREE R VALUE    /){
			print CHOTHIAPDB "$line\n";
		}elsif($identifier =~ "HETATM" && $line =~ /HOH/){
			#print CHOTHIAPDB "$line\n";
		}

	}
	close(PDBFILE);
}

sub readchothialight{
	@chothialight=[];
	my $var1=0;
	while(<CHOTHIALIGHT>){
		chomp;
		$var1++;
		$chothialight[$var1]=$_;
	}
}

sub readchothiaheavy{
	@chothiaheavy=[];
	my $var2=0;
	while(<CHOTHIAHEAVY>){
		chomp;
		$var2++;
		$chothiaheavy[$var2]=$_;
	}
}

sub removespaces {
	my $argument = $_[0] ;
	return ("") unless defined($argument);
	$argument =~ s/\s+$// ;
	$argument =~ s/^\s+// ;

	return $argument ;
}


