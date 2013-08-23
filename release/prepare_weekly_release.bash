#!/bin/bash
#The purpose of this script is to prepare a weekly release.  It WILL make git commits (no pushes), but I wouldn't suggest ever running it.

#It runs from the Rosetta folder (above main)

#check folder
for subdir in main tools demos
if [ ! -d $subdir ]
    echo "not running inside the Rosetta toplevel install directory; $subdir not found"
    exit 1
fi

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
