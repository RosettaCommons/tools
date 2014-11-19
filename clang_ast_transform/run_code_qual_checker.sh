#!/bin/bash

ME=$PWD
ROOT=/local/luki
export PYTHONPATH=$ROOT/tools/python_cc_reader

cd $ROOT/main-copy/source/src
python $ME/run_on_all_files_w_fork.py -e $ME/code_qual_check_file.sh -n 20
