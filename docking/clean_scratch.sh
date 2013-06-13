#!/bin/bash

#nodes=$(condor_status -long | grep ra | grep Machine |grep -v Client | cut -c 14-16 | uniq)
nodes=$(condor_status | grep @ | cut -d@ -f 2 | cut -d. -f 1 |sort |uniq)

echo nodes found: $nodes
echo user: $USER # should be found in the environment 

for node in $nodes
do
    echo Node $node
    if [ -d /scr/$node/$USER ]
    then
      cd /scr/$node/$USER/
      if [ "$1" = "erase" ]
      then
        rm -rf $2*
      fi
      du
    else
      echo Directory /scr/$node/$USER does not exist
    fi
done

cd /scr

echo usage: $(basename $0) [erase [subdir]] 