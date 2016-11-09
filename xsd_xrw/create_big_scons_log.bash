#!/usr/bin/env bash

set -x #echo commands
set -e #crash if any subcommand crashes

echo "usage: create_big_scons_log.bash #procs_to_compile_at"
rm -rf build/*
./scons.py -j$1  > big_scons_log
echo "big_scons_log should have your log file"
