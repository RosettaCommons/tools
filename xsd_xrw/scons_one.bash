#!/usr/bin/env bash

set -x #echo commands
set -e #crash if any subcommand crashes

if [ "$#" -ne 1 ]; then
    echo "usage: scons_one.bash file_to_compile;  partial path is also acceptable"
    exit
fi

echo "going to try this command:"
grep big_scons_log $1
grep big_scons_log $1 | source /dev/stdin
