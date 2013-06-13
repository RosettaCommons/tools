#!/usr/bin/perl

$label  = shift || die("usage: $0 column_label [datafile(s)]\n");

$first_line=<>; 


$column = index($first_line,$label);
#printf STDERR $first_line;
print "$label column =  $column\n";
