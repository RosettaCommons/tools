#!/bin/bash

ME=$PWD
ROOT=/data/rosetta
export PYTHONPATH=$ROOT/tools/python_cc_reader

cd $ROOT/main-copy/
#git reset --hard
#git stash apply

cd $ROOT/main/source/src
python $ME/run_on_all_files_w_fork.py -e $ME/rewrite_file.sh -n 8
