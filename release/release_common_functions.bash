#!/bin/bash
#THIS IS NOT A SCRIPT!  This is a collection of shared functions used in the release system.
#author: Steven Lewis, smlewi@gmail.com
#intended to be run only on Contador, but can be safely edited for whatever other machine is used to create the weekly Rosetta release.

#globally fail if any subcommand fails
set -e

#functions to read command line flags for scripts - this is adapted from http://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
echo "script takes these flags:\n
-d [no argument] debug mode\n
-j [integer argument] use this many processors MAXIMUM - may use less if machine load is high"

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:
debug=false
CONTADOR_MAX=24

while getopts "dj:" opt; do
    case "$opt" in
	d)  
	    debug=true
            ;;
	j)  
	    CONTADOR_MAX=$OPTARG
            ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

echo "debug=$debug, jobs max='$CONTADOR_MAX', Leftovers: $@"

#what is the name of this release?  Set up shared variables
week=$(date +%V)
year=$(date +%Y)
#Sergey's script to generate revID
#revID=sergeys_script.whatever
revID=revID

branch_name=$year/_$week/_revID
echo $branch_name " is branch name"

release_folder=Rosetta_$year.$week.revID
echo $release_folder " is release name"

#function call to "clean" a Rosetta install - removes all temp files, compiled files, etc
function simple_clean {

    if [ "$debug" = true ];
    then
	return 0 #do not perform clean in debug mode
    fi

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

#for "expected system load" when calculating how many processors to use
JOBS=0
function guess_load {
    uptime #for the user
    current_load=`uptime | awk -F '[ ,]' '{print $(NF-4)}'` #this parses "uptime" to grab the recent load
    floor_load=${current_load/.*} #I don't know what this does, but it floors the load value
    JOBS=$(($CONTADOR_MAX-(1+$floor_load))) #attempt number of jobs minus load ceiling
    echo "load was $current_load, attempting $JOBS"
} 

#determine our location, and bail out if this is the wrong place to be!
function check_folder {
    for subdir in main tools demos
    do
	if [ ! -d $subdir ]
	then
	    echo "not running inside the Rosetta toplevel install directory; $subdir not found"
	    exit 1
	fi
    done
}
