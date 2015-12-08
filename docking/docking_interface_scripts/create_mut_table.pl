#!/usr/bin/perl

#Author: Arvind Sivasubramanian
#Email: arvind_sivasub@yahoo.com

if( @ ARGV < 1 ){
die "Exiting. Not enough arguments...
Usage: $0 AXXB.C DYYYE.F ...
Script for creating mutants file for interface mode
A,B,D,E : single letter aa codes ( wt: A and D, mutations: B and E )
XX YYY: 2 or 3 digit residue number
C F : chain ids of the residues\n
Example: mutate.pl R21A.A S101A.B\n" ;
}

$nmutants = $#ARGV;
$printmutants = $nmutants + 1 ; #Number of mutations in the output

print "START\n";
print "$printmutants\n" ;

for ( $i = 0 ; $i <= $nmutants ; $i++ ) {

    $input = shift @ARGV ;
    @array = split (/\./ , $input ) ;

    $mutation = $array[0] ;
    $chain[$i] = $array[1] ;

    $length[$i] = $lngth = length ($mutation) ;
    $wt[$i] = substr ( $mutation , 0 , 1 ) ;
    $mutant[$i] = substr ( $mutation , $lngth-1 , 1 ) ;
    $resno[$i] = substr ( $mutation , 1 , $lngth - 2 ) ;
    print "1\n" ;
    printf "         %-4s  %s  %s  %s\n",$resno[$i],$chain[$i],$wt[$i],$mutant[$i];  

}
