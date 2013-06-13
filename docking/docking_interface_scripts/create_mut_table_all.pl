#!/usr/bin/perl

#Author Arvind Sivasubramanian
#Email: arvind_sivasub@yahoo.com

if( @ARGV < 2 ) {

die ( "Exiting. Not enough arguments...

Usage: $0 <interfacefile.int> <mutation: 1|2>
Script creates input file for the Rosetta Interface mode where each wt interface residue is mutated to most/all amino acids.
mutation = 1 : Mutate each interface residue to any one of the other 19 amino acids
mutation = 2 : Mutate to amino acids except to CYS and GLY
(Use findcontacts.py script to create .int file)
IMPORTANT: The actual number of mutations is printed as the last line of the output on the screen
Edit the input file for the interface/ddG mode to use this number after the \'START\' line\n" )
}


$interfacefile =  shift @ARGV ;
$mutation = shift @ARGV ;
$narguments = $#ARGV ;

%threetoone=("ALA",A,"CYS",C,"ASP",D,"GLU",E,"PHE",F,"GLY",G,"HIS",H,"ILE",I,"LYS",K,"LEU",L,"MET",M,"ASN",N,"PRO",P,"GLN",Q,"ARG",R,"SER",S,"THR",T,"VAL",V,"TRP",W,"TYR",Y);

@allaminoacids=("A","C","D","E","F", "G" , "H","I","K","L","M","N","P","Q","R","S","T","V","W","Y");
#No mutations to CYS or GLY
@allaminoacidsminus2=("A","D","E","F","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y");



@AAtomutate = @allaminoacids if ( $mutation == 1 ) ;
@AAtomutate = @allaminoacidsminus2 if ( $mutation == 2 ) ;


open ( interfacefile , $interfacefile ) || die ( "File $interfacefile not found\n" ) ;

$line = <interfacefile> ;
chop ( $line ) ;

until ( $line eq "" ) {

		unless ( $line =~ /Partner/ or $line =~ /GLY/ ) {
				$nmutants++ ;
				@array = split ( / +/ , $line ) ;

				$residuename[$nmutants] = $array[0] ;
				$chainid[$nmutants] = $array[1] ;
				$residuenumber[$nmutants] = $array[2] ;
		}

		$line = <interfacefile> ;
		chop ( $line ) ;
}

$number = @AAtomutate ;
$nmutations = $nmutants*$number ;
print "START\n" ;

print "$nmutations\n";

for($i=1;$i<=$nmutants;$i++) {

    $threelettercode=$residuename[$i] ;
    $wildtype=$threetoone{$threelettercode};

    foreach $mutant(@AAtomutate) {

				if($mutant ne $wildtype) {

						$number1++ ;
						print "1\n" ;
						printf "         %-4s  %s  %s  %s\n",$residuenumber[$i],$chainid[$i],$wildtype,$mutant ;
				}
		}

}
print "$number1\n" ;
