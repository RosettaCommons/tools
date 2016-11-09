#!/usr/bin/env bash

echo "usage: scons_one.bash file_to_compile;  partial path is also acceptable"
echo "going to try this command:"
grep big_scons_log $1
grep big_scons_log $1 | source /dev/stdin
