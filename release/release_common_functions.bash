#!/bin/bash
#THIS IS NOT A SCRIPT!  This is a collection of shared functions used in the release system.
#author: Steven Lewis, smlewi@gmail.com
#intended to be run only on Contador, but can be safely edited for whatever other machine is used to create the weekly Rosetta release.

#globally fail if any subcommand fails
set -e

#functions to read command line flags for scripts - this is adapted from http://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
echo "script takes these flags:
-d [no argument] debug mode
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
	    if (("$CONTADOR_MAX" > 24));
	    then
		echo "jobs variable (-j) CONTADOR_MAX cannot be greater than 24, it is $CONTADOR_MAX"
		exit 1
	    fi
            ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

echo "debug=$debug, jobs max='$CONTADOR_MAX', Leftovers: $@"

#what is the name of this release?  Set up shared variables
week=$(date +%V)
year=$(date +%Y)
#Sergey's script to generate revID needs to run post-git-happening
revID=revID

branch_name=$year\_$week\_revID
echo $branch_name " is tentative branch name"

release_folder=Rosetta_$year.$week.revID
echo $release_folder " is tentative release name"

#this function overwrites release_folder and branch_name with the real tag, GIVEN that you cd into main/ ahead of time.
function set_release_name {
    url_stem="http://benchmark.graylab.jhu.edu/data/get_commit_id/master/"
    git_hash=$(git rev-parse HEAD)
    echo "git hash for revID is $git_hash"
    wget -O revID_file $url_stem$git_hash
    revID=$(cat revID_file)
    rm revID_file

    branch_name=$year\_$week\_$revID
    echo $branch_name " is branch name"

    release_folder=Rosetta_$year.$week.$revID
    echo $release_folder " is release name"
}

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
    for subdir in main main/tools main/demos
    do
	if [ ! -d $subdir ]
	then
	    echo "not running inside the Rosetta toplevel install directory; $subdir not found"
	    exit 1
	fi
    done
}
