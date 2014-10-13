#!/bin/bash

ME=$PWD
ROOT=/local/luki
export PYTHONPATH=$ROOT/tools/python_cc_reader

cd $ROOT/main/source/src
python $ME/run_on_all_files_w_fork.py -e $ME/rewrite_file.sh -n 18

