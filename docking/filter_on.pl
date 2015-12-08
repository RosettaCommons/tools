#!/usr/bin/perl

#  if (@ARGV < 2) {
#      print STDERR "usage: $0 score lt/gt cutoff datafile\n";
#      exit;
#  }

$label  = shift;
$ltgt   = shift;
$cutoff = shift || die("usage: $0 column_label lt/gt cutoff datafile(s)\n");

if    ($ltgt eq "lt") { $lt++; $sym="<"; }
elsif ($ltgt eq "gt") {        $sym=">"; }
else { die("usagex: $0 column_label lt/gt cutoff datafile(s)\n"); }

$first_line=<>; 
$column = find_column($label,$first_line);
print $first_line;
printf STDERR "filtering on col $column = $label $sym $cutoff...";

$keepers=0;
$total_lines=0;
while (<>) {

    $total_lines++;

    @fields = split " ";
    next if (  ( $lt && $fields[$column] >= $cutoff) 
	     ||(!$lt && $fields[$column] <= $cutoff));

    $keepers++;
    print $_;

}

printf STDERR "kept $keepers out of $total_lines\n";
#print "\n";


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
    
