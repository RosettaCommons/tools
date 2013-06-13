#!/usr/bin/perl
## JJG revised 11/27/01

# now cuts off the worst 50%

$datafile=$ARGV[0];

# find the bump column 
$bump_label="fa_rep";
$bump_column=-1;
$header=`head -1 $datafile`;
@fields=split(/ +/,$header);
for ($i=0; $i<@fields; $i++) {
    if ($fields[$i] eq $bump_label) { $bump_column = $i; }
}
print STDERR "bump column $bump_column\n";
if ($bump_column == -1) { die("$bump_label column not found"); }

# find the 50% cutoff
$cutoff=`find_percent.pl $bump_column 50 $datafile`;
chomp($cutoff);
# used to bee minimum+100
#$cutoff=$minimum+100; # was *1.2, but not keeping very many in some cases; we really only want to get rid of egregious errors

print `filter_column.pl $bump_column $cutoff $datafile`;

