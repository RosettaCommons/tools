#!/usr/bin/perl
# jjg 3/14/02

if (@ARGV < 1) {
    print STDERR "usage: $0 column percent datafile(s)\n";
    print STDERR "finds the score at a given percentile within a column\n";
    exit;
}

$column  = shift;
$percent = shift;

while (<>) {

    @fields = split " ";

    if (/score/) {
	# header line
	printf STDERR "searching for $percent percent of $fields[$column]...";
#	push(@buf,$_);
	next;
    }

#    print $fields[0]," ",$fields[5],"\n";

    push(@values, $fields[$column]);
#    push(@buf,$_);

}

@sorted = sort { $a <=> $b } @values;
#print @sorted,"\n";

$index = $percent / 100.0 * scalar(@sorted);
#print $index,"\n";

$result = $sorted[$index];

printf STDERR "$result\n";
printf "$result\n";
