#!/bin/bash
#The purpose of this script is to delete integration tests from release packages.  It deletes tests that rely on unreleased code.  DO NOT RUN THIS SCRIPT EVER (except if you are preparing a release).  Run it from inside tests/integration.
#author: Steven Lewis, smlewi@gmail.com
#intended to be run only on Contador, but can be safely edited for whatever other machine is used to create the weekly Rosetta release.

#globally fail if any subcommand fails
set -e
source ./../../../tools/release/release_common_functions.bash

for directory in `ls tests`
do

#route around this one - it does bad things to our monster grep below (exposed by set -e; otherwise it just gets deleted and it's fine)
if ["$directory" = "app_exception_handling"];
    then
    continue
fi

#this monster finds all of the lines containing executable names from the tests' command files (first grep), strips the line down to only the part where the name is plus the flanking bin and binext string replacement junk (egrep), then also removes the string replacement junk (sed commands), then sorts and uniqs the list
for executable in `grep "%(bin)s" tests/$directory/command | grep -v MY_MINI_PROGRAM | egrep -o '%\(bin\)s/.*%\(binext\)s' | sed 's/%(bin)s\///g' | sed 's/.%(binext)s//g' | sort | uniq`
do
#echo $directory $executable

    if grep -q $executable ../../source/src/apps.src.settings
	then
	#if the executable is in the released build system, it's probably okay
	echo $executable "found in apps.src.settings" $directory "may be ok"
	else
	#if the executable is not found in the released build system, delete the test (as it won't run anyway)
	echo $executable "not found, REMOVE" $directory
	git rm -r tests/$directory #--dry-run
	if [ "$debug" = false ];
	    then
	    rm -r ref/$directory
	    else
	    echo "DEBUG MODE ACTIVATED: skipping filesystem ref deletion of autoremoved integration tests"
	fi
    fi
done

done
