#!/usr/bin/perl

# utility to take the output from charlie's histogram.pl and compare
# several files side-by-size
# jjgray 6/8/1

$nfiles=5;      # number of files
$labelsize=17;  # characters to save at the left of the output, from first file
$columnsize=15; # characters to take after that for data, from each file

@filedata = <>;
$filedatasize = @filedata;
$fsize = $filedatasize / $nfiles;

#print ("$fsize\n");

chop (@filedata);

### Header
print ("                 ");
for ($i=0;$i<$nfiles;$i++)
{
    print ("-----",substr(@filedata[$i*$fsize+0],6,4),"-----  ");
}
printf ("\n");

### Histograms
for ($line=1;$line<$fsize;$line++)
{
    print (substr(@filedata[0*$fsize+$line],1,$labelsize));
    for ($i=0;$i<$nfiles;$i++)
    {
	print (substr(@filedata[$i*$fsize+$line],$labelsize+1,$columnsize)," ");
    }
    printf ("\n");
}
