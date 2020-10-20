#!/bin/bash
# The purpose of this script is to convert standard rosetta checkout to its release form by removing unrelease parts of the code
# and deleting git information
# authosr: Steven Lewis, smlewi@gmail.com
#          Sergey Lyskov

#globally fail if any subcommand fails
set -e

source ./main/tools/release/release_common_functions.bash

check_folder #ensures we are in the right directory

ROSETTA=`pwd`

#update all repos
# if [ "$debug" = false ];
# then

#     #tools is a little weird to update, since this script updates overtop itself - so this really should be skipped while debugging
#     cd $ROSETTA/tools
#     pwd
#     git checkout master
#     git pull

#     cd $ROSETTA/demos
#     pwd
#     git checkout master
#     git pull

#     cd $ROSETTA/main
#     pwd
#     #notice this is release, not master, on the assumption that Sergey's machinery did some particular tests on each release push
#     git checkout release
#     git pull

#     cd $ROSETTA/documentation
#     pwd
#     git checkout master
#     git pull

# fi

#Make sure we have the latest of all the submodules
if [ "$debug" = false ];
then
    # cd $ROSETTA/tools
    # git submodule update --init --recursive

    # cd $ROSETTA/demos
    # git submodule update --init --recursive

    # cd $ROSETTA/documentation
    # git submodule update --init --recursive

    cd $ROSETTA/main
    # We have enough submodules in main that we may need to retry.
    # (On some systems Git/GitHub apparently throttles things.)
    until git submodule update --init --recursive; do
        echo ">>>>>>>>> Main submodule update did not complete successfully. Retrying. <<<<<<<<<"
        sleep 10
    done

fi


#prepare documentation
cd $ROSETTA/main/documentation
pwd

#git checkout -b weekly_releases/$branch_name #branch_name defined in release_common_functions
python remove_internal.py
#git commit -am "weekly release: remove 'internal' documentation"


#prepare tools
if [ "$debug" = false ];
then
    #tools is a little weird to update, since this script updates overtop itself - so this really should be skipped while debugging
    cd $ROSETTA/main/tools
    pwd
    #git checkout -b weekly_releases/$branch_name #branch_name defined in release_common_functions
fi

#prepare demos
cd $ROSETTA/main/demos
pwd
#git checkout -b weekly_releases/$branch_name #branch_name defined in release_common_functions

#prepare main with fresh & clean compile (for later itest references)

cd $ROSETTA
pwd

simple_clean #function call to simple_clean above
cd $ROSETTA/main/source
pwd

# guess_load
# #compile, run fixbb once (to get dunbrack binaries), delete fixbb, make a itest ref
# if [ "$debug" = false ];
# then
#     scons.py -j$JOBS bin mode=release

#     cd ../tests/integration
#     pwd
#     integration.py fixbb && rm -rf ref/ #get dunbrack binaries created so subsequent itest will have it premade
#     guess_load
#     integration.py -j $JOBS

# else
#     echo "DEBUG MODE ACTIVATED: skipping 'master' compile & test ref generation"
# fi

#we allow all these git commands to run in debug mode because it's on a branch anyway

#select release name
cd $ROSETTA/main
pwd
#set_release_name

cd $ROSETTA/main
pwd
#git checkout -b weekly_releases/$branch_name #branch_name defined in release_common_functions

for do_not_release in `cat $ROSETTA/main/tools/release/DONOTRELEASE.files`
do
    #sed -i -e '/DONOTRELEASE_TOP/,/DONOTRELEASE_BOTTOM/{/DONOTRELEASE_TOP/!{/DONOTRELEASE_BOTTOM/!d;}}' $do_not_release
    # <- does not work on Mac OS X due to outdated sed verion
    sed -i -e '/DONOTRELEASE_TOP/,/DONOTRELEASE_BOTTOM/s/.*/\/\/ DELETED BY RELEASE SCRIPT/' $do_not_release
done
#git commit -am "weekly release: stripping DONOTRELEASE-wrapped code"

cd $ROSETTA/main/source/src/devel
ls | grep -vE "init|svn_v" | xargs git rm -r
#git commit -m "weekly release: remove devel code"

cp init.release.cc init.cc
cp init.release.hh init.hh
#git commit -am "weekly release: devel/init.release"

cd $ROSETTA/main/source/src
pwd
cp devel.src.settings.release devel.src.settings
cp pilot_apps.src.settings.release pilot_apps.src.settings.all
#git commit -am "weekly release: overwrite *.src.settings with release versions"


# removing ui files
#cd $ROSETTA/main/source/ide && rm update_ui_project.py || true
cd $ROSETTA/main/source/src && git rm -rf ui || true


# Remove all of the pilot directories.
cd $ROSETTA/main/
find . -name pilot -type d -exec rm -r {} +
cd $ROSETTA/main/demos/
find . -name pilot -type d -exec rm -r {} +
#documentation and tools don't have pilot dirs, skip to speed things up.

# unneeded since Tim moved the documentation out; scheduled to be replaced with code that generates the static html wikimanual
# cd $ROSETTA/main/source/doc
# pwd
# git rm -r devel/ apps/pilot/
# git commit -m "weekly release: removing devel/pilot app documentation"

cd $ROSETTA/main/source/test/
pwd
cp devel.test.settings.release devel.test.settings
#git commit -am "weekly release: overwrite devel.test.settings with release version"

cd $ROSETTA/main/source/test/devel
ls | grep -vE "devel.cxxtest.hh" | xargs git rm -r
#git commit -m "weekly release: removing unit test devel"

# Remove Werror, etc from build settings
cd $ROSETTA/main/source/tools/build/
sed 's/^\(.*REMOVE FOR RELEASE\)/#\1/' basic.settings > basic.settings.release
mv -f basic.settings.release basic.settings

cd $ROSETTA/main/source/cmake/build/
sed 's/^\(.*REMOVE FOR RELEASE\)/#\1/' build.settings.cmake > basic.settings.cmake.release
mv -f basic.settings.cmake.release basic.settings.cmake

# Make the default build mode be release
cd $ROSETTA/main/source/tools/build/
sed 's/mode    = "debug"/mode    = "release"/' basic.options > basic.options.release
mv -f basic.options.release basic.options

# Enable path info scraping for released builds
cd $ROSETTA/main/source/tools/build/
cp site.settings.release site.settings


cd $ROSETTA/main/tests/integration/
pwd
#git rm -r tests/loop_creation tests/inverse_rotamer_remodel tests/splice_in tests/splice_out
rm -rf tests/loop_creation tests/inverse_rotamer_remodel

# if [ "$debug" = false ];
# then #must not do this rm - it is a filesystem rm of git ignored files
#     rm -r ref/loop_creation ref/inverse_rotamer_remodel ref/splice_in ref/splice_out
# else
#     echo "DEBUG MODE ACTIVATED: skipping filesystem ref deletion of manually removed integration tests"
# fi

#git commit -m "removing known-to-need-devel integration tests"
source $ROSETTA/main/tools/release/detect_itest_exes.bash
#git commit -m "deleting autoremoved integration tests"

#check compile
cd $ROSETTA/main/source/
# guess_load

# if [ "$debug" = false ];
# then
#     scons.py -j$JOBS bin mode=release

# else
#     echo "DEBUG MODE ACTIVATED: skipping post-releasification compile"
# fi

# echo "OK, the git branch should be ready...make sure that ^^ scons command worked, and look at the git history, then push with:"
# echo "cd main; git status; git log"
# echo "git push -u origin weekly_releases/$branch_name"


# Removing .git dirs -- We need to also remove the .git dirs of any submodule, hence the find/recursive

cd $ROSETTA
# Not just git directories, but the .git files for submodules, too.
find . -name '.git' -prune -exec rm -rf '{}' '+'
# We've collapsed everything, no need for submodules (which mess up the release git repo).
find . -name '.gitmodules' -type f -delete

echo "Conversion to release is done!"

exit
