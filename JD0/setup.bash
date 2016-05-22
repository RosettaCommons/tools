#!/bin/bash

cwd=`pwd`
mpth=$cwd/subdirs
echo $mpth

#split this file into 10 files, don't break lines, use numeric suffixes padded to 2 digits
split -a 2 -n l/10 --numeric-suffixes=1 rosetta_commands.list allpdbs_slice_

rm -r $mpth
mkdir $mpth

for subdir in `seq -w 10`
do
    echo $subdir
    mkdir $mpth/$subdir
    mv allpdbs_slice_$subdir $mpth/$subdir
    chmod 744 $mpth/$subdir/allpdbs_slice_$subdir
done

