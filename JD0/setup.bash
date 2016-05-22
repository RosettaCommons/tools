#!/bin/bash

cwd=`pwd`
mpth=$cwd/subdirs
echo $mpth

#split this file into 5 files, don't break lines, use numeric suffixes padded to 1 digit
#split -a 1 -n l/5 --numeric-suffixes=1 rosetta_commands.list allpdbs_slice_

#split this file into 10 files, don't break lines, use numeric suffixes padded to 2 digits
split -a 2 -n l/10 --numeric-suffixes=1 rosetta_commands.list allpdbs_slice_

#split this file into 100 files, don't break lines, use numeric suffixes padded to 3 digits
#split -a 3 -n l/100 --numeric-suffixes=1 rosetta_commands.list allpdbs_slice_

rm -r $mpth
mkdir $mpth

#for subdir in `seq -w 5`
for subdir in `seq -w 10`
#for subdir in `seq -w 100`

do
    echo $subdir
    mkdir $mpth/$subdir
    mv allpdbs_slice_$subdir $mpth/$subdir
    chmod 744 $mpth/$subdir/allpdbs_slice_$subdir
done

