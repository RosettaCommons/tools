#!/bin/bash
#The purpose of this script is to run in a cron job on the release-preparing machine, and set up the "release" branch for appropriate tests on Sergey's testing server.  In particular, it merges master to release around midnight on Thursdays.prepare a weekly release.  It WILL make a git ff-merge followed by a push, but should fail for non-ff merges requiring commits.
#It runs from the Rosetta folder (above main)
#author: Steven Lewis, smlewi@gmail.com

#globally fail if any subcommand fails
set -e

cd /media/scratch/smlewis/release_rosetta/Rosetta/main

git checkout master
git status

ideal_master_status="# On branch master
nothing to commit (working directory clean)"
real_master_status=$(git status)
if [ "$ideal_master_status" != "$real_master_status" ]
then
    pwd
    echo "update_release_branch_for_testing: abort, master was not clean"
    exit
fi

git pull

git checkout release
ideal_release_status="# On branch release
nothing to commit (working directory clean)"
real_release_status=$(git status)
if [ "$ideal_release_status" != "$real_release_status" ]
then
    pwd
    echo "update_release_branch_for_testing: abort, release was not clean"
    exit
fi

git merge --ff-only master

git push
