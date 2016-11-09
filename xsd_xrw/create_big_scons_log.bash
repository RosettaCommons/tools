#!/usr/bin/env bash

set -x #echo commands
set -e #crash if any subcommand crashes

if [ "$#" -ne 1 ]; then
    echo "usage: create_big_scons_log.bash #procs_to_compile_at"
    exit
fi

rm -rf build/*
#the grep will remove the scons startup stuff; those option file builds will conflict with some use cases
./scons.py -j$1  | grep -A 100000000 "scons: Building targets ..." > big_scons_log
echo "big_scons_log should have your log file"
