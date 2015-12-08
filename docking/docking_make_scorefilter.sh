# this script is called by rosetta.  it compiles the current
# scorefiles to determine where a certain percentage of scores lies,
# in order to filter out all but the top fraction of decoys
#
# jjg 12/15/01

percent=$1
homedir=$PWD
scdir=$(grep -E '^score' paths.txt | gawk {'print $2'})

cd $scdir
tmpf=tmp.sort
head -1 $(ls *sc |head -1) > $tmpf
cat *sc | grep -v filename | grep -v nat | grep -v inp | \
    sort -n +1 >> $tmpf

structs=$(wc $tmpf | gawk {'print $1'})
n=$[ $structs * $percent / 100 ]
nscore=$(head -$n $tmpf | tail -1 | gawk {'print $2'})

printf "$structs structures built\n"
printf "%10.1f  = score at $percent percent\n" $nscore

# contact filtered
filter_on.pl contact gt 1.9 $tmpf > tmp.sort.filtered
tmpf=tmp.sort.filtered

structs=$(wc $tmpf | gawk {'print $1'})
n=$[ $structs * $percent / 100 ]
nscore=$(head -$n $tmpf | tail -1 | gawk {'print $2'})

printf "$structs structures built with contact > 1.9\n"
printf "%10.1f  = score at $percent percent\n" $nscore
