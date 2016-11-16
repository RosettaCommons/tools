#!/usr/bin/env bash

set -x #echo commands
set -e #crash if any subcommand crashes

echo "git status suggests these:"
git status | grep "modified:" | awk '{print $2}'

echo "going to try these commands:"
git status | grep "modified:" | awk '{print $2}' | xargs -L 1 -I % grep % big_scons_log
git status | grep "modified:" | awk '{print $2}' | xargs -L 1 -I % grep % big_scons_log > temp_scons_stupid_mac_workaround
source temp_scons_stupid_mac_workaround
rm temp_scons_stupid_mac_workaround
