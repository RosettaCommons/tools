#!/bin/bash

# sorts a scorefile by a given score
# JJG 5/21/2004

if [ -z "$2" ]
then
    echo usage: $0 score scorefile [extra_flags]
    echo sorts a scorefile according to chosen score
    echo preserves the header line and automatically finds the right column
    exit;
fi

label=$1
scorefile=$2
shift;shift

col=$(findColumn.pl $label $scorefile | cut -f5 -d' ')

echo column: $col >2

head -1 $scorefile
sort -n +$col $@ $scorefile | grep -v filename
