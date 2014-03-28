ME=$PWD
ROOT=/data/rosetta
SRC=$ROOT/main/source/src

cd $SRC
export PYTHONPATH=$ROOT/tools/python_cc_reader

python $ME/run_on_all_files_w_fork.py -e $ME/rewrite_all.sh -n 8
