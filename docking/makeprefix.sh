#!/bin/bash

# converts a number to a rosetta prefix (0=aa, 1=ab... ZZ)

# if we already have a letter code, just return it
if [ $(echo $1 | grep -e "[[:alpha:]]") ]
then
    echo $1
    exit
fi


i=$[ $1 ]

letters=abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ

d1=$[ $i / 52 + 1 ]
d2=$[ $i - ( $d1 - 1 ) * 52 + 1 ]
#echo $d1 $d2

echo | gawk -v letters=$letters -v d1=$d1 -v d2=$d2 \
    '{print substr(letters,d1,1) substr(letters,d2,1); } '


# test code:

#  i=0
#  while [ $i -le 100 ]; do
#      echo $i  $(makeprefix.sh $i)
#      i=$[ i + 1 ]
#  done
