#!/usr/bin/perl

$label  = shift || die("usage: $0 column_label datafile(s)\n");

$first_line=<>; 
$column = find_column($label,$first_line);
#printf STDERR $first_line;
print "$label column =  $column\n";

sub find_column {
# find the column with specified label
    my ($label,$header) = @_;
    my $column = -1;
    my @fields=split(/ +/,$header);
    for ($i=0; $i<@fields; $i++) {
	if ($fields[$i] eq $label) { $column = $i; }
    }
    if ($column == -1) { die("$label column not found\nheader: $header"); }
    return $column;
}
    
