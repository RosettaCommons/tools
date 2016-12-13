#!/usr/bin/env bash

set -x #echo commands
set -e #crash if any subcommand crashes

if [ "$#" -ne 1 ]; then
    echo "usage: scons_one.bash file_to_compile;  partial path is also acceptable"
    exit
fi

echo "going to try this command:"
grep $1 big_scons_log
grep $1 big_scons_log  > temp_scons_stupid_mac_workaround
source temp_scons_stupid_mac_workaround
rm temp_scons_stupid_mac_workaround
