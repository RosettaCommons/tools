#!/bin/bash

cwd=`pwd`
mpth=$cwd/subdirs
echo $mpth

for subdir in `seq -w 10`
do
    echo $subdir
    cd $mpth/$subdir
    nohup ./allpdbs_slice_$subdir &

done

