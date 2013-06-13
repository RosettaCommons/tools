#!/usr/bin/perl

if ( $#ARGV < 2 ) {

    print "Not enough arguments. Exiting...\n" ;
    die ( "Usage : clustersbyscore.pl  <rmsfile> <clusteradius> <scorefile> \n " ) ;

}


$high = 1000000 ;
$rmsfile = shift @ARGV ;
$radius = shift @ARGV ;
$clusterfile = "tmp.clusters" ;
$scorefile = shift @ARGV ;
$outputfile = "centers" ;

print "Reading scorefile $scorefile ..." ;

open ( scores , $scorefile ) ;

$line = <scores> ;
chop ( $line ) ;

while ( $line ne '' ) {

( $pdbfile_w_leadingzeros , $score ) = split ( / +/ , $line ) ;

$pdbfile_wo_leadingzeros = $pdbfile_w_leadingzeros ;
$pdbfile_wo_leadingzeros =~ s/^0+// ;
$pdbfile_w_leadingzeros { $pdbfile_wo_leadingzeros } = $pdbfile_w_leadingzeros ;

$score { $pdbfile_wo_leadingzeros } = $score ;

$line = <scores> ;
chop ( $line ) ;


}

print "Done\n" ;

print "Clustering with $radius Angstrom radius...\n" ;

system " cat $rmsfile |rms2avglink1.csh $radius " ;

print "Processing clusters...\n" ;

open ( clusters , $clusterfile ) ;

$line = <clusters> ;
chop ( $line ) ;

while ( $line ne '' ) {

( $pdbfile , $clusternumber ) = split ( / +/ , $line ) ;

$members [ $clusternumber ] ++ ;

unless ( $seen [ $clusternumber ] ) {

    push @clusters , $clusternumber ;
    $best_score [ $clusternumber ] = $high ;
    $seen [ $clusternumber ] ++ ;



}

$filter_score = $score { $pdbfile } < $best_score [ $clusternumber ] ;

if ( $filter_score ) {

    $best_score [ $clusternumber ] = $score { $pdbfile } ;
    $center [ $clusternumber ] = $pdbfile ;

}


$line = <clusters> ;
chop ( $line ) ;


}


print "Printing output: centers centers.bysize centers.byscore ..." ;

open (outputfile , ">$outputfile" ) || die ( "File $outputfile cannot be found\n" ) ;

foreach $clusternumber ( @clusters ) {

    $center = $center [$clusternumber] ;

$best_score = $best_score [ $clusternumber ] ;

$pdbfile_w_leadingzeros = $pdbfile_w_leadingzeros { $center } ;


printf  outputfile "%-5d %5s %5d %8.2f\n" , $clusternumber , $pdbfile_w_leadingzeros , $members[ $clusternumber ] , $best_score ;


}


system "sort -nr +2 $outputfile > centers.bysize " ;
system "sort -n +3 $outputfile > centers.byscore " ; 

print "Done\n";
