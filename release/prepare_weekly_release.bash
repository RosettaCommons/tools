#!/bin/bash
#The purpose of this script is to prepare a weekly release.  It WILL make git commits (no pushes), but I wouldn't suggest ever running it.

#It runs from the Rosetta folder (above main)

#globally fail if any subcommand fails
set -e

#function call to "clean" a Rosetta install - removes all temp files, compiled files, etc
function simple_clean {
    if [ ! -d main ]
        then
        echo "simple_clean not running inside the Rosetta toplevel install directory; main not found"
        exit 1
    fi
    echo rm -r main/source/bin/*
    echo rm -r main/source/build/*
    echo rm main/source/.sconsign.dblite
    echo rm main/database/rotamer/bbdep02.May.sortlib.Dunbrack02.lib.bin
    echo rm main/database/rotamer/ExtendedOpt1-5/Dunbrack10.lib.bin
    echo rm main/source/.unit_test_results.yaml
    echo rm main/source/tools/build/user.options
    echo rm main/source/tools/build/user.settings
    echo rm -r main/tests/integration/new
    echo rm -r main/tests/integration/ref
    echo rm -r main/tests/integration/runtime_diffs.txt
    echo find . -name "*~" -exec rm {} \;
    echo find . -name "#*" -exec rm {} \;
    echo find . -name "*pyc" -exec rm {} \;
}

#!!!!!!!!!!!!!!!!!!!!!!!!!!!global variable !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#for "expected system load" when calculating how many processors to use
CONTADOR_MAX=24
JOBS=0
function guess_load {
    uptime #for the user
    current_load=`uptime | awk -F '[ ,]' '{print $(NF-4)}'` #this parses "uptime" to grab the recent load
    floor_load=${current_load/.*} #I don't know what this does, but it floors the load value
    JOBS=$(($CONTADOR_MAX-(1+$floor_load))) #attempt number of jobs minus load ceiling
    echo "load was $current_load, attempting $JOBS"
} 

#check folder
for subdir in main tools demos
do
    if [ ! -d $subdir ]
	then
	echo "not running inside the Rosetta toplevel install directory; $subdir not found"
	exit 1
    fi
done

ROSETTA=`pwd`

#update all repos
cd $ROSETTA/tools
pwd
echo git checkout master
echo git pull

cd $ROSETTA/demos
pwd
echo git checkout master
echo git pull

cd $ROSETTA/main
pwd
echo git checkout master
echo git pull

#prepare main with fresh & clean compile (for later itest references)
cd $ROSETTA
pwd
simple_clean #function call to simple_clean above
cd $ROSETTA/main/source
pwd
guess_load
#compile, run fixbb once (to get dunbrack binaries), delete fixbb, make a itest ref
echo scons.py -j$JOBS bin mode=release
cd ../tests/integration 
pwd
echo integration.py fixbb && rm -rf ref/
guess_load
echo integration.py -j $JOBS







exit




# for directory in `ls tests`
# do

# #this monster finds all of the lines containing executable names from the tests' command files (first grep), strips the line down to only the part where the name is plus the flanking bin and binext string replacement junk (egrep), then also removes the string replacement junk (sed commands), then sorts and uniqs the list
# for executable in `grep "%(bin)s" tests/$directory/command | grep -v MY_MINI_PROGRAM | egrep -o '%\(bin\)s/.*%\(binext\)s' | sed 's/%(bin)s\///g' | sed 's/.%(binext)s//g' | sort | uniq`
# do
# #echo $directory $executable

#     if grep -q $executable ../../source/src/apps.src.settings
# 	then
# 	#if the executable is in the released build system, it's probably okay
# 	echo $executable "found in apps.src.settings" $directory "may be ok"
# 	else
# 	#if the executable is not found in the released build system, delete the test (as it won't run anyway)
# 	echo $executable "not found, REMOVE" $directory
# 	git rm -r tests/$directory #--dry-run
# 	rm -r ref/$directory
#     fi
# done

# done
