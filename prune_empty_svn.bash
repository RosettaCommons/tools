#!/bin/bash

# This script finds empty folders in an svn repo. This is usually complicated because even "empty" folders contain a
# ".svn" directory and are not truly empty. If you uncomment some of the lines bellow it will print the full listing
# of each directory so you can confirm that all it contains is a .svn dir or just svn remove them for you.

# For those not fluent in bash/find/grep basically it iterates over each path that contains a directory called ".svn".
# For each one of those it does an "ls -a" but remove anything that has ".svn" in it. If the lenght of that string is
# zero (thats what the -z is for) it prints the name.

# It's a little round-about and kludgey but it works.
 
for i in `find ./ -type d -print -name ".svn" | grep -v "/\.svn"`
do
    if [ -z "`ls -A $i | grep -v '\.svn' `" ]
    then
        echo $i
        #ls -A $i
        #svn remove $i
    fi
done

