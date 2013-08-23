#!/bin/bash
#The purpose of this script is to prepare a weekly release.  It WILL make git commits (no pushes), but I wouldn't suggest ever running it.

#It runs from the Rosetta folder (above main)

#globally fail if any subcommand fails
set -e

#!!!!!!!!!!!!!!!!!!!!!!!!!what week is it?!!!!!!!!!!!!!!!!!!!!global
week=$(date +%V)


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

cd $ROSETTA/main
pwd
echo git checkout -b weekly_releases/2013-wk$week

cd $ROSETTA/main/source/src/devel
echo 'ls | grep -vE "init|svn_v" | xargs git rm -r'
echo git commit -m "weekly release: remove devel code"

echo cp init.release.cc init.cc
echo cp init.release.hh init.hh
echo git commit -am "weekly release: devel/init.release"

cd $ROSETTA/main/source/src
pwd
echo cp devel.src.settings.release devel.src.settings
echo cp pilot_apps.src.settings.release pilot_apps.src.settings.all
echo git commit -am "weekly release: overwrite *.src.settings with release versions"

cd $ROSETTA/main/source/src/apps
pwd
echo git rm -r curated pilot
echo git commit -m "weekly release: removing pilot apps"

cd $ROSETTA/main/source/doc
pwd
echo git rm -r devel/ apps/pilot/
echo git commit -m "weekly release: removing devel/pilot app documentation"

cd $ROSETTA/main/source/test/
pwd
echo cp devel.test.settings.release devel.test.settings
echo git commit -am "weekly release: overwrite devel.test.settings with release version"

cd $ROSETTA/main/source/test/devel
echo 'ls | grep -vE "devel.cxxtest.hh" | xargs git rm -r'
echo git commit -m "weekly release: removing unit test devel"

cd $ROSETTA/main/tests/integration/
pwd
echo git rm -r tests/loop_creation tests/inverse_rotamer_remodel tests/splice_in tests/splice_out 
echo rm -r ref/loop_creation ref/inverse_rotamer_remodel ref/splice_in ref/splice_out
echo git commit -m "removing known-to-need-devel integration tests"
echo $ROSETTA/tools/release/detect_itest_exes.bash
echo git commit -m "deleting autoremoved integration tests"

#THIS NEEDS MANUAL INTERVENTION
echo emacs -nw $ROSETTA/main/source/src/apps/public/design/mpi_msd.cc
echo git commit -am "weekly release: fixing the devel pilot app"

#check compile
cd $ROSETTA/main/source/
guess_load
echo scons.py -j$JOBS bin mode=release

echo "OK, the git branch should be ready...make sure that ^^ scons command worked, and look at the git history, then push with:"
echo "????"

exit
