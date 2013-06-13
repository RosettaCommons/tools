#!/usr/bin/perl

#Author: Arvind Sivasubramanian
#Email: arvind_sivasub@yahoo.com

if ( @ARGV < 1 ) {

die ( "Exiting. Not enough arguments...
Usage: $0 <energycutoff>

For each decoy (pdbname.pdb), this script parses the corresponding Rosetta Interface mode 
output file (pdbname.out) and outputs the number of '\Hits\'. A '\Hit\' results when there
is agreement between the experimentally determined effect of a mutation and the Rosetta
Interface prediction.

Empirical \'Hit\' definitions ( Kortemme and Baker, PNAS 2002)

Rosetta ddG < 0 : Binding increase
0 < Rosetta ddG < energycutoff : Neutral
Rosetta ddG > energycutoff : Binding loss

Edit the script to add entries to the following three lists

\@bindingloss: Mutations that are experimentally known to result in binding loss
\@neutral: Mutations that are experimentally known not to affect binding one way or another
\@mutation_subset: Subset of the bindingloss/neutral mutations for which the analysis is required

For instance \@bindingloss = ( N682A , L685A ) and so on
\n" ) ;

}
$energycutoff = shift @ARGV;
$energycutoff = 1.0 if ($energycutoff eq '');

@bindingloss = (N682A,L685A,Y688A) ; 
@neutral = (K684A) ;
@mutation_subset= (@bindingloss,@neutral) ;


foreach $mutation(@bindingloss){
    $result{$mutation}="Loss";
}

foreach $mutation(@neutral){
    $result{$mutation}="Neutral";
}

@files=<*.out>;
$filecount=0;

foreach $file(@files){

    $filecount++ ;
    $fileid="mutfile".$filecount;

    open(fileid,$file);
    $line=<fileid>;

    $mutantno=0;
    $nhits=$nmiss=0;
    while($line ne "") {
				$mutantno++ ;
				chop($line) ;
				@contents=split(/ +/,$line) ;

				if($contents[0] eq "DDG_BIND") {

						$energy=$contents[13] ;
						$resno_chain=$contents[23] ;
						$wt_mut=$contents[24] ;

						$resno = substr ( $resno_chain , 0 , length ( $resno_chain ) - 1 ) ;
						$wt = substr ( $wt_mut , 0 , 1 ) ;
						$mut = substr ( $wt_mut , 2 , 1 ) ;


						$mutation=$wt.$resno.$mut ;

						foreach $mutant(@mutation_subset) {

								if($mutation eq $mutant) {

										if ($energy > $max_energy{$mutation}) {
												$max_energy{$mutation}=$energy  ;
												$max_decoy{$mutation}=$file ;
										}

										if ( $energy < 0 ) {
												$decoyresult = "Negative" ;
										} elsif ( $energy > 0 && $energy < $energycutoff ) {
												$decoyresult = "Neutral" ;
										} else {
												$decoyresult = "Loss" ;
										}

										if($decoyresult eq $result{$mutation}) {
												$finalresult="Hit" ;
												$nhits++ ;
										}
										else {
												$finalresult="Miss" ;
												$nmiss++ ;
										}

										$numberofhits{$mutation}++ if( $finalresult eq "Hit") ;
								}
						}
				}

				$line=<fileid> ;
		}
    $nmutations=$nhits+$nmiss ;
    $frachits=$nhits*1.0/$nmutations if($nmutations != 0);
    $fracmiss=1-$frachits;
    $hitmissratio=$nhits*1.0/$nmiss if($nmiss != 0);
    $pdbfile = substr($file,0,length($file)-4).".pdb";
    printf "%19s     %4.2f   %4.2f   %3s %3s %3s\n",$pdbfile,$frachits,$fracmiss,$nhits,$nmiss,$nmutations ;
}

