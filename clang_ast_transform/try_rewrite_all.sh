#!/bin/bash

ME=$PWD
ROOT=/data/rosetta
export PYTHONPATH=$ROOT/tools/python_cc_reader

cd $ROOT/main-copy/
git reset --hard

rsync -a $ROOT/main/source/src/ $ROOT/main-copy/source/src/

cd $ROOT/main/source/src
python $ME/run_on_all_files_w_fork.py -e $ME/rewrite_file.sh -n 8

cd $ROOT/main-copy/
patch -p0 --merge < ../main-patches/all.patch
patch -p0 --merge < ../main-patches/settings.patch

#for P in $ROOT/main-patches/*.patch; do
#	patch -p1 --merge < $P
#done

cd $ROOT/main-copy/source
scons cxx=clang mode=debug -k -j 4
