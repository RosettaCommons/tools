#!/usr/bin/env bash

set -x #echo commands
set -e #crash if any subcommand crashes

if [ "$#" -ne 1 ]; then
    echo "usage: create_big_scons_log.bash #procs_to_compile_at"
    exit
fi

rm -rf build/*
./scons.py -j$1  > big_scons_log
echo "big_scons_log should have your log file"
