#!/bin/bash

cwd=`pwd`
mpth=$cwd/subdirs
echo $mpth

#for subdir in `seq -w 5`
for subdir in `seq -w 10`
#for subdir in `seq -w 100`

do
    echo $subdir
    cd $mpth/$subdir
    nohup ./allpdbs_slice_$subdir &

done

