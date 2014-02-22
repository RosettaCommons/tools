#!/bin/bash
#The purpose of this script is to move and test the weekly release.
#It runs from the Rosetta folder (above main), and contains HARDCODED paths to the release-prep machine

#globally fail if any subcommand fails
set -e

#!!!!!!!!!!!!!!!!!!!!!!!!!what week is it?!!!!!!!!!!!!!!!!!!!!global
week=$(date +%V)
year=$(date +%Y)

#function call to "clean" a Rosetta install - removes all temp files, compiled files, etc
function simple_clean {
    if [ ! -d main ]
        then
        echo "simple_clean not running inside the Rosetta toplevel install directory; main not found"
        exit 1
    fi
    set +e #globally, we need exit-on-error, but it's ok if these rms fail to find targets
    rm -r main/source/bin/*
    rm -r main/source/build/*
    rm main/source/.sconsign.dblite
    rm main/database/rotamer/bbdep02.May.sortlib.Dunbrack02.lib.bin
    rm main/database/rotamer/ExtendedOpt1-5/Dunbrack10.lib.bin
    rm main/database/rotamer/bbdep02.May.sortlib-correct.12.2010.Dunbrack02.lib.bin
    rm main/source/.unit_test_results.yaml
    rm main/source/tools/build/user.options
    rm main/source/tools/build/user.settings
    rm -r main/tests/integration/new
    rm -r main/tests/integration/ref
    rm -r main/tests/integration/runtime_diffs.txt
    find . -name "*~" -exec rm {} \;
    find . -name "#*" -exec rm {} \;
    find . -name "*pyc" -exec rm {} \;
    set -e #return to exit-on-error
}

function deep_clean {
simple_clean

rm -rf demos/.git
rm -rf tools/.git
rm -rf main/.git
rm main/.gitmodules
find . -name ".gitignore" -exec rm {} \;

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

cd -P  .. #this is above Rosetta/
pwd
#tar Rosetta
tar -cf Rosetta_unstripped_release.tar Rosetta/

cd -P  /media/scratch/smlewis/release_holding_area
pwd
mv /media/scratch/smlewis/release_rosetta/Rosetta_unstripped_release.tar .

#..untar a copy
tar -xf Rosetta_unstripped_release.tar
release_folder=Rosetta_wk$week\_$year
mv Rosetta/ $release_folder
rm Rosetta_unstripped_release.tar
cd -P  $release_folder
pwd

#clean the copy
deep_clean

#tar up the candidate release
cd -P ..
pwd
tar -czf $release_folder.tar.gz $release_folder
cd -P $release_folder
pwd

guess_load
#compile, run fixbb once (to get dunbrack binaries), delete fixbb, make a itest ref
cd -P  main/source/
pwd
scons.py -j$JOBS bin mode=release
scons.py -j$JOBS
scons.py -j$JOBS cat=test
test/run.py -j$JOBS

cd -P  ../tests/integration
pwd
integration.py fixbb && rm -rf ref/
integration.py -j $JOBS

mv ref/ new/
cp -ar /media/scratch/smlewis/release_rosetta/Rosetta/main/tests/integration/ref/ .

echo "check the unit and itests, if it compiled it's probably good"
echo "cd `pwd` && integration.py --compareonly --fulldiff"

exit
