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
    rm -r main/source/bin/*
    rm -r main/source/build/*
    rm main/source/.sconsign.dblite
    rm main/database/rotamer/bbdep02.May.sortlib.Dunbrack02.lib.bin
    rm main/database/rotamer/ExtendedOpt1-5/Dunbrack10.lib.bin
    rm main/source/.unit_test_results.yaml
    rm main/source/tools/build/user.options
    rm main/source/tools/build/user.settings
    rm -r main/tests/integration/new
    rm -r main/tests/integration/ref
    rm -r main/tests/integration/runtime_diffs.txt
    find . -name "*~" -exec rm {} \;
    find . -name "#*" -exec rm {} \;
    find . -name "*pyc" -exec rm {} \;
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
git checkout master
git pull

cd $ROSETTA/demos
pwd
git checkout master
git pull

cd $ROSETTA/main
pwd
git checkout master
git pull

#prepare main with fresh & clean compile (for later itest references)
cd $ROSETTA
pwd
simple_clean #function call to simple_clean above
cd $ROSETTA/main/source
pwd
guess_load
#compile, run fixbb once (to get dunbrack binaries), delete fixbb, make a itest ref
scons.py -j$JOBS bin mode=release
cd ../tests/integration 
pwd
integration.py fixbb && rm -rf ref/
guess_load
integration.py -j $JOBS

cd $ROSETTA/main
pwd
git checkout -b weekly_releases/2013-wk$week

cd $ROSETTA/main/source/src/devel
ls | grep -vE "init|svn_v" | xargs git rm -r
git commit -m "weekly release: remove devel code"

cp init.release.cc init.cc
cp init.release.hh init.hh
git commit -am "weekly release: devel/init.release"

cd $ROSETTA/main/source/src
pwd
cp devel.src.settings.release devel.src.settings
cp pilot_apps.src.settings.release pilot_apps.src.settings.all
git commit -am "weekly release: overwrite *.src.settings with release versions"

cd $ROSETTA/main/source/src/apps
pwd
git rm -r curated pilot
git commit -m "weekly release: removing pilot apps"

cd $ROSETTA/main/source/doc
pwd
git rm -r devel/ apps/pilot/
git commit -m "weekly release: removing devel/pilot app documentation"

cd $ROSETTA/main/source/test/
pwd
cp devel.test.settings.release devel.test.settings
git commit -am "weekly release: overwrite devel.test.settings with release version"

cd $ROSETTA/main/source/test/devel
ls | grep -vE "devel.cxxtest.hh" | xargs git rm -r
git commit -m "weekly release: removing unit test devel"

cd $ROSETTA/main/tests/integration/
pwd
git rm -r tests/loop_creation tests/inverse_rotamer_remodel tests/splice_in tests/splice_out 
rm -r ref/loop_creation ref/inverse_rotamer_remodel ref/splice_in ref/splice_out
git commit -m "removing known-to-need-devel integration tests"
$ROSETTA/tools/release/detect_itest_exes.bash
git commit -m "deleting autoremoved integration tests"

#THIS NEEDS MANUAL INTERVENTION
emacs -nw $ROSETTA/main/source/src/apps/public/design/mpi_msd.cc
git commit -am "weekly release: fixing the devel pilot app"

#check compile
cd $ROSETTA/main/source/
guess_load
scons.py -j$JOBS bin mode=release

echo "OK, the git branch should be ready...make sure that ^^ scons command worked, and look at the git history, then push with:"
echo "????"

exit
