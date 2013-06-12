#!/usr/bin/perl

if ( @ARGV < 4 ) {

print "Not enough arguments. Exiting...
Usage: $0 <chainID> <start_res> <end_res> <pdbfile>\n" ;

die ( ) ;

}

$chnid=shift @ARGV;
$lowres=shift @ARGV;
$highres= shift @ARGV;
$pdbfilename= shift @ARGV;

#Testing for insertion code in residue number
$lowres_insertion = '';
$highres_insertion = ' ';
$last_char = substr( $lowres, -1, 1);
if( $last_char =~ m/[A-Za-z]/ ) {
		$lowres_insertion = $last_char;
}
$last_char = substr( $highres, -1, 1);
if( $last_char =~ m/[A-Za-z]/ ) {
		$highres_insertion = $last_char;
}
$highres_insertion_flag = 1;

open(pdbfile,$pdbfilename);

$line=<pdbfile>;
chop($line);
while($line ne ""){

    if ($line =~ /ATOM/){
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
				$bfactor = substr($line,60,6);
				
				$residueno =~ s/^\s+//;
				$residueno =~ s/\s+$//;
				
				if($chain eq $chnid){
						if($residueno >= $lowres){
								if( ($lowres_insertion eq '') || ($insert_code eq $lowres_insertion)){
										$lowres_insertion = '';
										if($residueno <= $highres){
												if( $highres_insertion_flag || ($highres_insertion eq $insert_code)) { 
														print "$line\n";
												}
												if( ($residueno == $highres) && ($highres_insertion eq $insert_code) ) {
														$highres_insertion_flag = 0;
												}
										}
								}
						}
				}
		}
		
    $line=<pdbfile>;
    chop($line);
}
