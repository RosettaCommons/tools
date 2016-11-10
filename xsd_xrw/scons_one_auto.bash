#!/usr/bin/env bash

set -x #echo commands
set -e #crash if any subcommand crashes

echo "git status suggests these:"
git status | grep "modified:" | awk '{print $2}'

echo "going to try these commands:"
git status | grep "modified:" | awk '{print $2}' | xargs -L 1 -I % grep % big_scons_log
git status | grep "modified:" | awk '{print $2}' | xargs -L 1 -I % grep % big_scons_log | source /dev/stdin
