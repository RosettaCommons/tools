#!/bin/bash
#The purpose of this script is to move and test the weekly release.
#It runs from the Rosetta folder (above main), and contains HARDCODED paths to the release-prep machine
#author: Steven Lewis, smlewi@gmail.com
#intended to be run only on Contador, but can be safely edited for whatever other machine is used to create the weekly Rosetta release.

#globally fail if any subcommand fails
set -e

source ./main/tools/release/release_common_functions.bash

function deep_clean {

    if [ "$debug" = true ];
    then
	return 0 #do not perform clean in debug mode
    fi

    simple_clean
    rm -rf demos/.git
    rm -rf tools/.git
    rm -rf main/.git
    rm main/.gitmodules
    find . -name ".gitignore" -exec rm {} \;

}

check_folder #ensures we are in the right directory

#set release name (must happen in git-land)
cd main
pwd
set_release_name
cd ..

cd -P  .. #this is above Rosetta/
pwd
#tar Rosetta
tar -cf Rosetta_unstripped_release.tar Rosetta/

cd -P  /media/scratch/smlewis/release_holding_area
pwd
mv /media/scratch/smlewis/release_rosetta/Rosetta_unstripped_release.tar .

#..untar a copy
tar -xf Rosetta_unstripped_release.tar

#release_folder set in release_common_functions.bash
mv Rosetta/ $release_folder
rm Rosetta_unstripped_release.tar
cd -P $release_folder
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
if [ "$debug" = false ];
then
    #compile
    scons.py -j$JOBS bin mode=release
    scons.py -j$JOBS
    scons.py -j$JOBS cat=test

    #unit test
    test/run.py -j$JOBS

    #integration test
    cd -P  ../tests/integration
    pwd
    integration.py fixbb && rm -rf ref/
    integration.py -j $JOBS

    #set up test diffs
    mv ref/ new/
    cp -ar /media/scratch/smlewis/release_rosetta/Rosetta/main/tests/integration/ref/ .
else
    echo "DEBUG MODE ACTIVATED: not actually recompiling/testing candidate release tarball"
fi

echo "check the unit and itests, if it compiled it's probably good"
echo "cd `pwd` && integration.py --compareonly --fulldiff"

exit
