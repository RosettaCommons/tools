#!/bin/bash

cpus=$(nproc)
rosetta_location=~/Rosetta/main

if [ ! -d $rosetta_location ]; then
    echo "Rosetta location '$rosetta_location' doesn't exist"
    exit
fi

cd $rosetta_location/source/cmake/build_debug
if [[ $# -gt 0 ]]; then
    if [[ $1 =~ .*?[.]test$ ]]; then
	nice -n 15 ninja -j $cpus $2
    else
	nice -n 15 ninja -j $cpus $2_symlink &> /tmp/compilation_build_output.txt &
    fi
else
    nice -n 15 ninja -j $cpus &> /tmp/compilation_build_output.txt &
fi

echo "Release build:"
cd $rosetta_location/source/cmake/build_release
if [[ $# -gt 1 ]]; then
    if [[ $2 =~ .*?[.]test$ ]]; then
	echo "Skipping release build of $2"
    else
	nice -n 5 ninja -j $cpus $2_symlink
    fi
else
    nice -n 5 ninja -j $cpus
fi
return_val=$?
wait

if [ $return_val = 0 ]; then
    echo "Debug build:"
    cat /tmp/compilation_build_output.txt
fi

rm /tmp/compilation_build_output.txt &> /dev/null
exit $return_val
